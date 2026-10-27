#!/usr/bin/env python3
"""Summarize raw or filtered DDA peptide identifications from idXML."""
import argparse
from collections import Counter, defaultdict
import pyopenms as oms
from QCMetrics import Metrics, finite, text


def label_of(hit):
    """Return the target/decoy class of a peptide hit."""
    if not hit.metaValueExists('target_decoy'):
        return 'unknown'
    value = text(hit.getMetaValue('target_decoy')).lower()
    return value if value in ('target', 'decoy') else 'unknown'


def score_rows(identifications, filtered):
    """Collect selected hit scores and attributes by target/decoy class."""
    scores = defaultdict(list)
    signed_errors = defaultdict(list)
    absolute_errors = defaultdict(list)
    charges = defaultdict(Counter)
    missing_charges = Counter()
    missing_errors = Counter()
    score_types, directions = set(), set()
    total_hits, represented_spectra = 0, set()
    sequences = set()
    for identification in identifications:
        hits = list(identification.getHits())
        total_hits += len(hits)
        if not hits:
            continue
        score_types.add(text(identification.getScoreType()))
        higher = identification.isHigherScoreBetter()
        directions.add(higher)
        spectrum_ref = (text(identification.getMetaValue('spectrum_reference'))
                        if identification.metaValueExists('spectrum_reference') else '')
        if spectrum_ref:
            merge_index = (text(identification.getMetaValue('id_merge_index'))
                           if identification.metaValueExists('id_merge_index') else '')
            represented_spectra.add((merge_index, spectrum_ref))
        valid = [(finite(hit.getScore()), hit) for hit in hits]
        valid = [(score, hit) for score, hit in valid if score is not None]
        selected = hits if filtered else ([
            (max(valid, key=lambda pair: pair[0]) if higher else
             min(valid, key=lambda pair: pair[0]))[1]] if valid else [])
        for hit in selected:
            label = label_of(hit)
            score = finite(hit.getScore())
            sequence = text(hit.getSequence().toString())
            if sequence:
                sequences.add(sequence)
            for group in ('all', label):
                scores[group].append(score)
            charge = hit.getCharge()
            if charge > 0:
                for group in ('all', label):
                    charges[group][charge] += 1
            else:
                for group in ('all', label):
                    missing_charges[group] += 1
            observed_mz = finite(identification.getMZ())
            theoretical_mz = (finite(hit.getSequence().getMZ(charge))
                              if charge > 0 else None)
            if observed_mz is not None and theoretical_mz is not None and theoretical_mz > 0:
                error = (observed_mz - theoretical_mz) / theoretical_mz * 1_000_000
                for group in ('all', label):
                    signed_errors[group].append(error)
                    absolute_errors[group].append(abs(error))
            else:
                for group in ('all', label):
                    missing_errors[group] += 1
    return (scores, signed_errors, absolute_errors, charges, missing_charges,
            missing_errors, score_types, directions, total_hits, represented_spectra,
            sequences)


def summarize(args):
    """Read idXML and write search or validation QC."""
    protein_ids = []
    peptide_ids = oms.PeptideIdentificationList()
    oms.IdXMLFile().load(args.input, protein_ids, peptide_ids)
    filtered = args.stage == 'filtered'
    layer = 'validated_identification' if filtered else 'raw_identification'
    metrics = Metrics(args.dataset_id, layer, args.input, sample_name=args.sample_name,
                      search_engine=args.search_engine)
    (scores, errors, absolute, charges, missing_charges, missing_errors,
     score_types, directions, total_hits, spectra, sequences) = score_rows(peptide_ids, filtered)
    metrics.add('peptide_identification_count', len(peptide_ids))
    metrics.add('accepted_psm_count' if filtered else 'candidate_psm_count', total_hits)
    metrics.add('distinct_peptide_count', len(sequences))
    if not filtered:
        metrics.add('top_hit_count', len(scores['all']))
        if spectra:
            metrics.add('spectrum_with_candidate_count', len(spectra))
        else:
            metrics.unavailable('spectrum_with_candidate_count')
        metrics.unavailable('searched_ms2_spectrum_count')
        metrics.unavailable('spectrum_match_coverage', 'fraction')
    if len(score_types) == 1:
        metrics.add('score_type', next(iter(score_types)), '')
    elif score_types:
        raise ValueError('Mixed score types in one idXML; separate them before QC')
    if len(directions) == 1:
        metrics.add('score_higher_is_better', int(next(iter(directions))), 'boolean')
    for group in ('all', 'target', 'decoy', 'unknown'):
        metrics.distribution('validated_score' if filtered else 'search_score', scores[group],
                             group=group)
        if not filtered:
            metrics.distribution('precursor_mass_error_signed', errors[group], 'ppm', group)
            metrics.distribution('precursor_mass_error_absolute', absolute[group], 'ppm', group)
            metrics.add('precursor_mass_error_missing_count', missing_errors[group], group=group)
            metrics.add('precursor_charge_missing_count', missing_charges[group], group=group)
            total = sum(charges[group].values())
            for charge, count in sorted(charges[group].items()):
                metrics.add('precursor_charge_count', count, group=f'{group}:z{charge}')
                metrics.add('precursor_charge_fraction', count / total, 'fraction',
                            f'{group}:z{charge}')
    metrics.write(args.output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('input', 'output', 'dataset-id', 'sample-name', 'search-engine'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--stage', choices=('search', 'filtered'), required=True)
    summarize(parser.parse_args())
