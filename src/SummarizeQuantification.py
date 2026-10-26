#!/usr/bin/env python3
"""Summarize OpenMS ProteinQuantifier peptide or protein abundance TSV."""
import argparse
import csv
import math
import re
from pathlib import Path
from QCMetrics import Metrics, finite


def table_lines(handle, comments):
    """Yield table lines while collecting OpenMS provenance comments."""
    for line in handle:
        if line.startswith('#'):
            comments.append(line)
        elif line.strip():
            yield line


def summarize(args):
    """Report feature coverage, missingness, and positive abundance distributions."""
    metrics = Metrics(args.dataset_id, 'quantification', args.input)
    with open(args.input, newline='', encoding='utf-8') as handle:
        comments = []
        reader = csv.DictReader(table_lines(handle, comments), delimiter='\t')
        columns = reader.fieldnames or []
        feature_column = 'protein' if args.level == 'protein' else 'peptide'
        if feature_column not in columns:
            raise ValueError(f'{args.input} lacks {feature_column} column')
        abundance_columns = [name for name in columns if name.startswith('abundance_')
                             or name.startswith('abundance|')]
        if not abundance_columns:
            raise ValueError(f'{args.input} has no ProteinQuantifier abundance columns')
        values = {name: [] for name in abundance_columns}
        features = {name: set() for name in abundance_columns}
        row_count = 0
        for row in reader:
            row_count += 1
            for name in abundance_columns:
                value = finite(row[name])
                values[name].append(value)
                if value is not None and value > 0:
                    features[name].add(row[feature_column])
    metrics.add(args.level + '_table_row_count', row_count)
    assay_names = {}
    for comment in comments:
        for number, run in re.findall(r"fgroup(\d+): '([^']+)'", comment):
            assay_names[number] = run
    for name in abundance_columns:
        assay = re.fullmatch(r'abundance_fgroup(\d+)_label\d+', name)
        if assay:
            metrics.context['sample_name'] = assay_names.get(assay.group(1), '')
        elif name.startswith('abundance|'):
            metrics.context['sample_name'] = Path(name.split('|')[1]).stem
        column = values[name]
        group = f'{args.level}:{name}'
        missing = sum(value is None for value in column)
        zeros = sum(value == 0 for value in column)
        metrics.add('quantified_feature_count', len(features[name]), group=group)
        metrics.add('abundance_missing_count', missing, group=group)
        metrics.add('abundance_zero_count', zeros, group=group)
        metrics.add('abundance_missing_fraction', missing / row_count if row_count else None,
                    'fraction', group, 'available' if row_count else 'unavailable')
        metrics.distribution('abundance_log10',
                             (math.log10(value) if value is not None and value > 0
                              else None for value in column), 'log10 intensity', group)
    metrics.write(args.output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ('input', 'output', 'dataset-id'):
        parser.add_argument('--' + argument, required=True)
    parser.add_argument('--level', choices=('protein', 'peptide'), required=True)
    summarize(parser.parse_args())
