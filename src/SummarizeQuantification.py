#!/usr/bin/env python3
"""Summarize ProteinQuantifier abundance by LFQ assay or TMT channel."""
import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path
from QCMetrics import Metrics, finite


def table_lines(handle, comments):
    """Yield table lines while collecting OpenMS provenance comments."""
    for line in handle:
        if line.startswith('#'):
            comments.append(line)
        elif line.strip():
            yield line


def read_acquisitions(path):
    """Map exact mzML and sample names to sample/fraction pairs."""
    sources = defaultdict(set)
    pairs = set()
    with open(path, newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle)
        if not {'sample_name', 'fraction_id', 'path'} <= set(reader.fieldnames or []):
            raise ValueError('mzML CSV needs sample_name,fraction_id,path')
        for row in reader:
            sample, fraction = row['sample_name'].strip(), row['fraction_id'].strip()
            source = Path(row['path'].strip())
            pair = (sample, fraction)
            if not sample or not fraction or not source.name or pair in pairs:
                raise ValueError(f'Invalid or duplicate mzML sample/fraction: {pair}')
            pairs.add(pair)
            for key in (sample, source.name, source.stem):
                sources[key].add(pair)
    if not pairs:
        raise ValueError('mzML CSV has no acquisitions')
    return sources, pairs


def read_channels(path, samples):
    """Read TMT descriptions without assuming reporter-channel order."""
    descriptions = defaultdict(set)
    seen = set()
    with open(path, newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle)
        if not {'plex_id', 'channel', 'description'} <= set(reader.fieldnames or []):
            raise ValueError('TMT channel CSV needs plex_id,channel,description')
        for row in reader:
            plex = row['plex_id'].strip()
            channel = row['channel'].strip()
            description = row['description'].strip()
            if plex not in samples or not channel or (plex, channel) in seen:
                raise ValueError(f'Invalid or duplicate TMT channel: {plex}/{channel}')
            seen.add((plex, channel))
            if description:
                descriptions[plex].add(description)
    if not seen:
        raise ValueError('TMT channel CSV has no channels')
    return descriptions


def read_assay_mapping(path, pairs):
    """Read explicit mappings for ambiguous abundance headers."""
    if path == 'none':
        return {}
    mapping = {}
    required = {'level', 'abundance_column', 'sample_name', 'fraction_id',
                'reporter_channel', 'group'}
    samples = {pair[0] for pair in pairs}
    with open(path, newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle)
        if not required <= set(reader.fieldnames or []):
            raise ValueError('Assay mapping CSV needs ' + ','.join(sorted(required)))
        for row in reader:
            key = (row['level'].strip(), row['abundance_column'].strip())
            sample, fraction = row['sample_name'].strip(), row['fraction_id'].strip()
            if (key in mapping or key[0] not in ('peptide', 'protein') or not key[1]
                    or sample not in samples
                    or (fraction and (sample, fraction) not in pairs)):
                raise ValueError(f'Invalid or duplicate assay mapping: {key}')
            mapping[key] = (sample, fraction, row['reporter_channel'].strip(),
                            row['group'].strip())
    return mapping


def abundance_header(name):
    """Decode a current ProteinQuantifier abundance header."""
    parts = name.split('|')
    if len(parts) not in (3, 4, 5) or parts[0] != 'abundance':
        raise ValueError(f'Unsupported abundance header: {name}')
    if not parts[1] or not re.fullmatch(r'ch[1-9][0-9]*', parts[2]):
        raise ValueError(f'Invalid abundance source/channel: {name}')
    if len(parts) == 5 and not re.fullmatch(r'fgroup[1-9][0-9]*', parts[4]):
        raise ValueError(f'Invalid abundance assay suffix: {name}')
    return parts[1], parts[2], parts[3] if len(parts) >= 4 else ''


def resolve_assay(name, level, mode, sources, pairs, channels, mapping):
    """Resolve sample, fraction, and TMT label without guessing."""
    source, reporter, description = abundance_header(name)
    if (level, name) in mapping:
        sample, fraction, mapped_reporter, group = mapping[level, name]
        if mode == 'TMT' and mapped_reporter != reporter:
            raise ValueError(f'Assay mapping has wrong reporter channel for {name}')
        if mode == 'LFQ' and (mapped_reporter or group not in ('', 'all')):
            raise ValueError(f'LFQ assay mapping has TMT dimensions for {name}')
    else:
        candidates = set()
        for key in (source, Path(source).name, Path(source).stem):
            candidates.update(sources.get(key, ()))
        if not candidates and mode == 'TMT' and description:
            matches = {sample for sample, labels in channels.items()
                       if description in labels}
            if len(matches) == 1:
                candidates = {pair for pair in pairs if pair[0] in matches}
        sample_names = {pair[0] for pair in candidates}
        if len(sample_names) != 1:
            raise ValueError(f'Cannot map {name} to one sample; supply assay_mapping_csv')
        sample = sample_names.pop()
        fraction = next(iter(candidates))[1] if len(candidates) == 1 else ''
        group = description if mode == 'TMT' else 'all'
    if mode == 'TMT':
        if group and group not in channels[sample]:
            raise ValueError(f'TMT sample label {group!r} is absent from channel CSV for {sample}')
        channel = reporter
    else:
        group, channel = 'all', ''
    if level == 'protein':
        fraction = ''
    if fraction and (sample, fraction) not in pairs:
        raise ValueError(f'Assay {name} refers to unknown fraction {sample}/{fraction}')
    return sample, fraction, channel, group


class QuantificationSummarizer:
    """Summarize each ProteinQuantifier abundance assay at peptide or group level."""

    def __init__(
        self,
        input_path,
        dataset_id,
        level,
        mode,
        mzml_csv,
        channel_csv,
        assay_mapping_csv,
    ):
        self.input_path = input_path
        self.level = level
        self.mode = mode
        self.mzml_csv = mzml_csv
        self.channel_csv = channel_csv
        self.assay_mapping_csv = assay_mapping_csv
        self.metrics = Metrics(
            dataset_id=dataset_id,
            layer='quantification',
            source=input_path,
        )

    def load_assay_context(self):
        """Read sample, fraction, channel, and explicit assay mappings."""
        sources, pairs = read_acquisitions(self.mzml_csv)
        samples = {sample for sample, _ in pairs}
        channels = (read_channels(self.channel_csv, samples)
                    if self.mode == 'TMT' else {})
        mapping = read_assay_mapping(self.assay_mapping_csv, pairs)
        return sources, pairs, channels, mapping

    def read_abundance_table(self):
        """Collect abundance values and distinct positive feature identities."""
        comments = []
        with open(self.input_path, newline='', encoding='utf-8') as handle:
            reader = csv.DictReader(table_lines(handle, comments), delimiter='\t')
            columns = reader.fieldnames or []
            if self.level not in columns:
                raise ValueError(f'{self.input_path} lacks {self.level} column')
            abundance_columns = [name for name in columns if name.startswith('abundance')]
            if not abundance_columns or len(abundance_columns) != len(set(abundance_columns)):
                raise ValueError(f'{self.input_path} has no unique abundance columns')
            for name in abundance_columns:
                abundance_header(name)
            values = {name: [] for name in abundance_columns}
            features = {name: set() for name in abundance_columns}
            for row in reader:
                for name in abundance_columns:
                    value = finite(row[name])
                    values[name].append(value)
                    if value is not None and value > 0 and row[self.level]:
                        features[name].add(row[self.level])
        return values, features, comments

    @staticmethod
    def reconcile_fraction_groups(sources, comments):
        """Connect OpenMS fraction-group aliases to manifest source names."""
        for comment in comments:
            for number, run in re.findall(r"fgroup(\d+): '([^']+)'", comment):
                for key in (run, Path(run).name, Path(run).stem):
                    sources[f'fgroup{number}'].update(sources.get(key, ()))

    def summarize_assays(self, values, features, sources, pairs, channels, mapping):
        """Write distinct identity counts and positive log10 distributions."""
        prefix = 'peptide' if self.level == 'peptide' else 'protein_group'
        for name in values:
            sample, fraction, _, group = resolve_assay(
                name=name,
                level=self.level,
                mode=self.mode,
                sources=sources,
                pairs=pairs,
                channels=channels,
                mapping=mapping,
            )
            self.metrics.context.update(
                sample_name=sample,
                fraction_id=fraction,
            )
            self.metrics.add(
                name=f'quantified_{prefix}_count',
                value=len(features[name]),
                group=group,
            )
            self.metrics.distribution(
                name=f'{prefix}_abundance_log10',
                values=(math.log10(value) if value is not None and value > 0 else None
                        for value in values[name]),
                unit='log10 intensity',
                group=group,
            )

    def run(self, output_path):
        """Resolve assays, summarize their abundance, and write QC rows."""
        sources, pairs, channels, mapping = self.load_assay_context()
        values, features, comments = self.read_abundance_table()
        self.reconcile_fraction_groups(sources, comments)
        self.summarize_assays(
            values=values,
            features=features,
            sources=sources,
            pairs=pairs,
            channels=channels,
            mapping=mapping,
        )
        self.metrics.write(output_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ('input', 'output', 'dataset-id', 'mzml-csv', 'channel-csv'):
        parser.add_argument('--' + argument, required=True)
    parser.add_argument('--level', choices=('protein', 'peptide'), required=True)
    parser.add_argument('--quantification-mode', choices=('TMT', 'LFQ'), required=True)
    parser.add_argument('--assay-mapping-csv', default='none')
    args = parser.parse_args()
    summarizer = QuantificationSummarizer(
        input_path=args.input,
        dataset_id=args.dataset_id,
        level=args.level,
        mode=args.quantification_mode,
        mzml_csv=args.mzml_csv,
        channel_csv=args.channel_csv,
        assay_mapping_csv=args.assay_mapping_csv,
    )
    summarizer.run(args.output)
