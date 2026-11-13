#!/usr/bin/env python3
"""Summarize cohort protein groups in OpenMS idXML or consensusXML."""
import argparse
import pyopenms as oms
from QCMetrics import Metrics, finite, text


class ProteinGroupsSummarizer:
    """Summarize one cohort protein-inference artifact."""

    def __init__(self, input_path, dataset_id, route):
        self.input_path = input_path
        self.route = route
        self.metrics = Metrics(
            dataset_id=dataset_id,
            layer='protein_inference',
            source=input_path,
        )
        self.groups = []
        self.score_types = set()

    def load_identifications(self):
        """Load ProteinProphet idXML or ProteinInferencer consensusXML."""
        protein_ids = []
        if self.route == 'ProteinProphet':
            peptide_ids = oms.PeptideIdentificationList()
            oms.IdXMLFile().load(self.input_path, protein_ids, peptide_ids)
        else:
            consensus = oms.ConsensusMap()
            oms.ConsensusXMLFile().load(self.input_path, consensus)
            protein_ids = consensus.getProteinIdentifications()
        return protein_ids

    def collect_groups(self, protein_ids):
        """Keep group members, probabilities, and input score types."""
        for identification in protein_ids:
            self.score_types.add(text(identification.getScoreType()))
            selected = (identification.getProteinGroups() if self.route == 'ProteinProphet'
                        else identification.getIndistinguishableProteins())
            for group in selected:
                members = list(group.accessions)
                if not members:
                    raise ValueError('Protein group has no members')
                self.groups.append((members, finite(group.probability)))

    def summarize_membership(self):
        """Count inferred groups and their member proteins."""
        self.metrics.add('protein_group_count', len(self.groups))
        self.metrics.distribution(
            name='proteins_per_group',
            values=(len(members) for members, _ in self.groups),
        )

    def summarize_scores(self):
        """Use probabilities for ProteinProphet and scores for ProteinInferencer."""
        if self.route == 'ProteinProphet':
            self.metrics.distribution(
                name='protein_group_probability',
                values=(value for _, value in self.groups),
                unit='probability',
            )
        else:
            if len(self.score_types) == 1:
                self.metrics.add(
                    name='protein_group_score_type',
                    value=next(iter(self.score_types)),
                    unit='',
                )
            self.metrics.distribution(
                name='protein_group_score',
                values=(value for _, value in self.groups),
                unit='score',
            )

    def run(self, output_path):
        """Read inference groups and write cohort QC metrics."""
        self.collect_groups(self.load_identifications())
        self.summarize_membership()
        self.summarize_scores()
        self.metrics.write(output_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('input', 'output', 'dataset-id'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--route', choices=('ProteinProphet', 'ProteinInferencer'), required=True)
    args = parser.parse_args()
    summarizer = ProteinGroupsSummarizer(
        input_path=args.input,
        dataset_id=args.dataset_id,
        route=args.route,
    )
    summarizer.run(args.output)
