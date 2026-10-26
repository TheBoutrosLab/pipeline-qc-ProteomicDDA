#!/usr/bin/env python3
"""Summarize one DDA mzML fraction with pyOpenMS."""
import argparse
import pyopenms as oms
from QCMetrics import Metrics, finite


def summarize(args):
    """Measure recorded duration, scan rates, cycles, and TIC."""
    experiment = oms.MSExperiment()
    oms.MzMLFile().load(args.input, experiment)
    metrics = Metrics(args.dataset_id, 'mzML', args.input, sample_name=args.sample_name,
                      fraction_id=args.fraction_id)
    times, ms1, ms2, tic = [], [], [], []
    cycle_counts = []
    current_ms2 = 0
    for spectrum in experiment.getSpectra():
        rt = finite(spectrum.getRT())
        if rt is None or (times and rt < times[-1]):
            raise ValueError('mzML has missing or nonmonotonic retention time')
        times.append(rt)
        if spectrum.getMSLevel() == 1:
            if ms1:
                cycle_counts.append(current_ms2)
            ms1.append(rt)
            current_ms2 = 0
        elif spectrum.getMSLevel() == 2:
            ms2.append(rt)
            current_ms2 += 1
        _, intensities = spectrum.get_peaks()
        tic.append(float(sum(intensities)) if len(intensities) else None)
    if not times:
        raise ValueError('mzML has no spectra')
    duration = times[-1] - times[0]
    metrics.add('acquisition_duration', duration / 60, 'min')
    for name, count in [('scan_count', len(times)), ('ms1_scan_count', len(ms1)),
                        ('ms2_scan_count', len(ms2))]:
        metrics.add(name, count)
        rate_name = name.replace('count', 'rate')
        if duration > 0:
            metrics.add(rate_name, count / duration, 'scan/s')
        else:
            metrics.unavailable(rate_name, 'scan/s')
    metrics.distribution('ms1_cycle_time', (b - a for a, b in zip(ms1, ms1[1:])), 's')
    metrics.distribution('ms2_scans_per_cycle', cycle_counts)
    metrics.distribution('total_ion_current', tic, 'intensity')
    midpoint = len(tic) // 2
    metrics.distribution('total_ion_current_early', tic[:midpoint], 'intensity')
    metrics.distribution('total_ion_current_late', tic[midpoint:], 'intensity')
    metrics.write(args.output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ('input', 'output', 'sample-name', 'fraction-id', 'dataset-id'):
        parser.add_argument('--' + argument, required=True)
    summarize(parser.parse_args())
