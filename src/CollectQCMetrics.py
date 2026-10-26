#!/usr/bin/env python3
"""Combine DDA QC metric TSVs into one deterministic table."""
import argparse
import csv
from QCMetrics import FIELDS


def collect(args):
    """Validate schema, combine, and sort metrics."""
    combined = []
    for path in args.inputs:
        with open(path, newline='', encoding='utf-8') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            if reader.fieldnames != FIELDS:
                raise ValueError(f'Unexpected QC metric columns in {path}: {reader.fieldnames}')
            combined.extend(reader)
    combined.sort(key=lambda row: tuple(row[key] for key in
                  ('dataset_id', 'layer', 'sample_name', 'fraction_id', 'search_engine',
                   'inference_route', 'metric_name', 'group', 'source_file')))
    with open(args.output, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, FIELDS, delimiter='\t')
        writer.writeheader()
        writer.writerows(combined)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', nargs='+', required=True)
    parser.add_argument('--output', required=True)
    collect(parser.parse_args())
