# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
[unreleased]: https://github.com/TheBoutrosLab/pipeline-qc-ProteomicDDA/compare/v1.1.0...HEAD
