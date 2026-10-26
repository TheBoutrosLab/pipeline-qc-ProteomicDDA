#!/usr/bin/env python3
"""Summarize cohort protein groups in OpenMS idXML or consensusXML."""
import argparse
import pyopenms as oms
from QCMetrics import Metrics, finite, text


def summarize(args):
    """Load the selected inference artifact and summarize group membership."""
    protein_ids = []
    if args.route == 'ProteinProphet':
        peptide_ids = oms.PeptideIdentificationList()
        oms.IdXMLFile().load(args.input, protein_ids, peptide_ids)
    else:
        consensus = oms.ConsensusMap()
        oms.ConsensusXMLFile().load(args.input, consensus)
        protein_ids = consensus.getProteinIdentifications()
    groups = []
    score_types = set()
    for identification in protein_ids:
        score_types.add(text(identification.getScoreType()))
        selected = (identification.getProteinGroups() if args.route == 'ProteinProphet'
                    else identification.getIndistinguishableProteins())
        for group in selected:
            members = list(group.accessions)
            if not members:
                raise ValueError('Protein group has no members')
            groups.append((members, finite(group.probability)))
    metrics = Metrics(args.dataset_id, 'protein_groups', args.input,
                      inference_route=args.route)
    metrics.add('protein_group_count', len(groups))
    metrics.distribution('proteins_per_group', (len(members) for members, _ in groups))
    if args.route == 'ProteinProphet':
        metrics.distribution('protein_group_probability',
                             (value for _, value in groups), 'probability')
    else:
        if len(score_types) == 1:
            metrics.add('protein_group_score_type', next(iter(score_types)), '')
        metrics.distribution('protein_group_score',
                             (value for _, value in groups), 'score')
    metrics.write(args.output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('input', 'output', 'dataset-id'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--route', choices=('ProteinProphet', 'ProteinInferencer'), required=True)
    summarize(parser.parse_args())
