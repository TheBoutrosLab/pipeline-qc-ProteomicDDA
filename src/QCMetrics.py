"""Shared numeric summaries and long-form DDA QC output."""
import csv
import math

FIELDS = ['dataset_id', 'acquisition_mode', 'layer', 'sample_name', 'fraction_id',
          'search_engine', 'inference_route', 'source_process', 'metric_name',
          'group', 'value', 'unit', 'status', 'source_file', 'calculation_version']
QUANTILES = (5, 25, 50, 75, 95)


def text(value):
    """Convert OpenMS byte strings to text."""
    return value.decode() if isinstance(value, bytes) else str(value)


def finite(value):
    """Return finite float or None."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def percentile(sorted_values, percent):
    """Calculate an interpolated percentile."""
    position = (len(sorted_values) - 1) * percent / 100
    left, right = math.floor(position), math.ceil(position)
    return sorted_values[left] + (sorted_values[right] - sorted_values[left]) * (position - left)


class Metrics:
    """Collect stage metrics with explicit grain and availability."""
    def __init__(self, dataset_id, layer, source, **dimensions):
        self.context = {'dataset_id': dataset_id, 'acquisition_mode': 'DDA',
                        'layer': layer, 'sample_name': '', 'fraction_id': '',
                        'search_engine': '', 'inference_route': '', 'source_process': '',
                        'source_file': str(source), 'calculation_version': '1.0.0'}
        self.context.update(dimensions)
        self.rows = []

    def add(self, name, value, unit='count', group='all', status='available'):
        """Add one metric, including explicit unavailable status."""
        self.rows.append(dict(self.context, metric_name=name, group=group,
                              value='' if value is None else value, unit=unit,
                              status=status))

    def unavailable(self, name, unit='', group='all'):
        """Record unavailable source data without treating it as zero."""
        self.add(name, None, unit, group, 'unavailable')

    def distribution(self, name, values, unit='', group='all'):
        """Write count, missing count, and compact percentiles."""
        values = list(values)
        numbers = sorted(number for value in values if (number := finite(value)) is not None)
        self.add(name + '_count', len(numbers), group=group)
        self.add(name + '_missing_count', len(values) - len(numbers), group=group)
        for q in QUANTILES:
            if numbers:
                self.add(f'{name}_p{q:02d}', percentile(numbers, q), unit, group)
            else:
                self.unavailable(f'{name}_p{q:02d}', unit, group)

    def write(self, path):
        """Write metric rows to a TSV."""
        with open(path, 'w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, FIELDS, delimiter='\t')
            writer.writeheader()
            writer.writerows(self.rows)
