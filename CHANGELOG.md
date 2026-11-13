# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.2.0] - 2026-11-13

### Added

- Add a schematic of the DDA QC metric and MS1 visualization workflows.

### Changed

- Match the DIA QC table's 12 columns and their order while standardizing shared metric names.
- Count raw top-hit peptides separately from all validated peptides; align score metadata and unknown target/decoy label counts.
- Summarize ProteinQuantifier proteins as quantified protein groups, with TMT biological labels in `group` and LFQ samples in `sample_name`.
- Refresh LFQ and TMT integration fixtures to verified metapipeline 1.2.7 and ProteinQuantifier 1.3.0 outputs.
- Organize QC readers, summarizers, matrix extraction, plotting, and metric collection into focused classes without changing their output.

### Removed

- Remove DDA-only QC table columns for acquisition mode, quantification mode, inference route, reporter channel, assay ID, and calculation version.
- Drop DDA early/late TIC, unavailable searched-spectrum denominator/coverage, cohort quantification row totals, and separate abundance zero/missing metrics.
- Drop compatibility with older ProteinQuantifier abundance headers and previous QC output names.

## [1.1.1] - 2026-10-28

### Added

- Count pre-validation MSBooster PSMs by missed cleavages, enzyme cleavage specificity, and precursor charge for all, target, and decoy groups.

## [1.1.0] - 2026-10-27

### Added

- Summarize published MSBooster PIN and PSMFeatureExtractor idXML features from an optional manifest, with separate target and decoy distributions.
- Record the producing feature process in each QC metric row.
- Refresh LFQ and TMT ProteinProphet integration inputs to published metapipeline 1.0.1 outputs, including PSMFeatureExtractor idXML files.

### Changed

- Use stable pipeline layer names and sort the merged QC table in pipeline order.

## [1.0.0] - 2026-10-26

### Added

- Initial Nextflow DSL2 DDA quality-control pipeline for completed metapipeline outputs, with versioned Boutros Lab configuration, container profiles, resource allocation, and a combined long-form QC metrics TSV.
- Per-fraction mzML acquisition summaries and offline interactive 3D MS1 peak maps generated through separate pyOpenMS extraction and Plotly rendering processes.
- Unfiltered search and final filtered identification summaries from idXML, including target and decoy score and mass-error distributions.
- Cohort protein-group summaries for ProteinProphet and ProteinInferencer, and peptide and protein abundance summaries from ProteinQuantifier.
- UCLA nf-test integration cases for LFQ and TMT ProteinProphet outputs, with case-specific metapipeline inputs.
- Run, input, output, and testing documentation; individual process logs remain available after workflow completion.

[1.0.0]: https://github.com/TheBoutrosLab/pipeline-qc-ProteomicDDA/releases/tag/v1.0.0
[1.1.0]: https://github.com/TheBoutrosLab/pipeline-qc-ProteomicDDA/compare/v1.0.0...v1.1.0
[1.1.1]: https://github.com/TheBoutrosLab/pipeline-qc-ProteomicDDA/compare/v1.1.0...v1.1.1
[1.2.0]: https://github.com/TheBoutrosLab/pipeline-qc-ProteomicDDA/compare/v1.1.1...v1.2.0
[unreleased]: https://github.com/TheBoutrosLab/pipeline-qc-ProteomicDDA/compare/v1.2.0...HEAD
