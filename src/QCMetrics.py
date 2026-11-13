"""Shared numeric summaries and long-form DDA QC output."""
import csv
import math

FIELDS = ['dataset_id', 'layer', 'source_process', 'sample_name', 'search_engine',
          'fraction_id', 'metric_name', 'group', 'value', 'unit', 'status', 'source_file']
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
        self.context = {'dataset_id': dataset_id, 'layer': layer, 'source_process': '',
                        'sample_name': '', 'search_engine': '', 'fraction_id': '',
                        'source_file': str(source)}
        self.context.update(dimensions)
        self.rows = []

    def add(self, name, value, unit='count', group='all', status='available'):
        """Add one metric, including explicit unavailable status."""
        row = {
            **self.context,
            'metric_name': name,
            'group': group,
            'value': '' if value is None else value,
            'unit': unit,
            'status': status,
        }
        self.rows.append(row)

    def unavailable(self, name, unit='', group='all'):
        """Record unavailable source data without treating it as zero."""
        self.add(
            name=name,
            value=None,
            unit=unit,
            group=group,
            status='unavailable',
        )

    def distribution(self, name, values, unit='', group='all'):
        """Write count, missing count, and compact percentiles."""
        values = list(values)
        numbers = sorted(number for value in values if (number := finite(value)) is not None)
        self.add(
            name=name + '_count',
            value=len(numbers),
            group=group,
        )
        self.add(
            name=name + '_missing_count',
            value=len(values) - len(numbers),
            group=group,
        )
        for q in QUANTILES:
            if numbers:
                self.add(
                    name=f'{name}_p{q:02d}',
                    value=percentile(numbers, q),
                    unit=unit,
                    group=group,
                )
            else:
                self.unavailable(
                    name=f'{name}_p{q:02d}',
                    unit=unit,
                    group=group,
                )

    def unavailable_distribution(self, name, unit='', group='all'):
        """Mark an entire distribution unavailable when its source field is absent."""
        for suffix in ('count', 'missing_count', 'p05', 'p25', 'p50', 'p75', 'p95'):
            self.unavailable(
                name=f'{name}_{suffix}',
                unit='' if suffix.endswith('count') else unit,
                group=group,
            )

    def write(self, path):
        """Write metric rows to a TSV."""
        with open(path, 'w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, FIELDS, delimiter='\t')
            writer.writeheader()
            writer.writerows(self.rows)
