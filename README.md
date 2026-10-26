# pipeline-qc-ProteomicDDA

This Nextflow DSL2 pipeline collects QC metrics from completed `metapipeline-ProteomicDDA` output. It reads each original mzML fraction, unfiltered search idXML for each sample and engine, the final filtered validation idXML for each sample or TMT plex, one cohort protein-inference artifact, and the ProteinQuantifier peptide and protein tables. It also creates an interactive 3D MS1 peak map for every mzML fraction.

## How to run

1. Initialize Git submodules with `git submodule update --init --recursive`.
2. Copy [config/template.config](config/template.config) to a run-specific config and set all input paths, `inference_route`, `dataset_id`, and `output_dir`.
3. Run with Nextflow 24.10.6 and the v1 syntax parser:

   ```bash
   export NXF_VER=24.10.6
   export NXF_SYNTAX_PARSER=v1
   nextflow run . -c /path/to/qc.config -profile apptainer
   ```

The `standard` and `docker` profiles use Docker; `apptainer` uses Apptainer. The OpenMS processes use `ghcr.io/theboutroslab/pyopenms:3.5.0`; the HTML plot uses `multiqc/multiqc:v1.35`, which provides Plotly. Set `work_dir` to scratch storage when processing large mzML files.

## Inputs

Three CSV manifests identify the same set of samples or TMT plexes. Paths in these CSVs may be absolute or relative to the CSV directory. Multiple mzML fractions and multiple search engines per sample are supported.

**mzML CSV** (`mzml_csv`):

```csv
sample_name,fraction_id,path
sample_01,1,/data/sample_01_fraction_1.mzML
sample_01,2,/data/sample_01_fraction_2.mzML
```

**Unfiltered search idXML CSV** (`search_idxml_csv`):

```csv
sample_name,search_engine,idXML
sample_01,MSFragger,/data/sample_01_MSFragger.idXML
```

**Final filtered idXML CSV** (`filtered_idxml_csv`):

```csv
sample_name,idXML
sample_01,/data/sample_01_MSFragger_filter.idXML
```

For TMT, `plex_id,idXML` is accepted in the filtered CSV. The search CSV can be the metapipeline's published `validation.csv`, and the filtered CSV can be its published `quant-idxml.csv`. Each mzML sample must occur in both identification manifests. Use one final consensus filtered idXML per sample when searches from multiple engines were merged; do not sum it with per-engine filtered counts.

| Parameter | Description |
| --- | --- |
| `dataset_id` | Cohort identifier used in the output path and metric filenames. |
| `mzml_csv` | Per-fraction mzML manifest. |
| `search_idxml_csv` | Unfiltered search idXML manifest, keyed by sample and search engine. |
| `filtered_idxml_csv` | Final filtered idXML manifest, keyed by sample or TMT plex. |
| `inference_route` | `ProteinProphet` or `ProteinInferencer`. |
| `protein_groups` | Cohort filtered ProteinProphet idXML or ProteinInferencer consensusXML. |
| `protein_quant_csv` | ProteinQuantifier protein abundance table. |
| `peptide_quant_csv` | ProteinQuantifier peptide abundance table. |
| `output_dir` | Root directory for QC results. |
| `work_dir` | Optional Nextflow work directory. |
| `pyopenms_version` | pyOpenMS image tag; default `3.5.0`. |
| `docker_image_plotly` | Plotly image; default `multiqc/multiqc:v1.35`. |
| `base_resource_update` | Optional process CPU and memory overrides. |

ProteinProphet grouping uses the filtered idXML passed to ProteinQuantifier. The protein quantitation child currently publishes it only when `save_intermediate_files = true`. ProteinInferencer grouping uses the inferred consensusXML; choose the protein-FDR-filtered version if that was passed downstream. The two ProteinQuantifier tables are tab-delimited despite their `.csv` extension. Their comment lines and `abundance_fgroup*_label*` or detailed `abundance|...|ch*` columns are parsed directly. Run names from the fraction-group comment populate `sample_name` when available; the abundance-column ID remains in `group` to preserve TMT label and assay identity.

## Pipeline steps

| Process | Input | Output |
| --- | --- | --- |
| `summarize_mzML` | One original mzML fraction | Acquisition duration, MS1/MS2 scan counts and rates, survey cycle timing, MS2 scans per cycle, and whole-run/early/late TIC summaries |
| `extract_MS1_matrix` | One mzML fraction | Bounded RT/m/z/intensity peak matrix and selection metadata |
| `plot_MS1_landscape` | Matrix and metadata | Offline interactive 3D MS1 stick spectrum |
| `summarize_search_idXML` | One raw search idXML per sample and engine | Candidate and top-hit counts, score, signed/absolute ppm mass error, target/decoy, and precursor charge summaries |
| `summarize_filtered_idXML` | One final filtered idXML per sample or plex | Accepted PSMs, distinct retained peptides, and validated score summaries |
| `summarize_ProteinGroups` | One cohort inference artifact | Group count and size, plus ProteinProphet probability or ProteinInferencer score summaries when present |
| `summarize_ProteinQuantifier` | Peptide and protein abundance tables | Quantified feature counts, zero/missing counts, missing fractions, and log10 positive intensity summaries per assay column |
| `collect_QCMetrics` | All stage TSVs | Combined `<dataset_id>_QC_metrics.tsv` |

A continuous distribution is saved as its nonmissing count, missing count, and 5th, 25th, 50th, 75th, and 95th percentiles. `group=all` includes target and decoy top hits; target and decoy also have separate summaries. Charge states have counts and fractions. The searchable-spectrum denominator is not present in the raw idXML when unmatched spectra are omitted, so `searched_ms2_spectrum_count` and `spectrum_match_coverage` are explicitly marked `unavailable`. A recorded mzML acquisition duration is not necessarily the programmed LC gradient duration.

The MS1 matrix includes up to 120 strongest positive peaks from each of 350 evenly spaced MS1 scans within m/z 300–1200. The HTML embeds Plotly and works offline. Its colors encode within-run log-intensity percentiles: yellow (0–25), orange (25–50), red (50–75), purple (75–95), and blue (95–100). Colors do not represent the same absolute intensity across runs.

## Outputs

Results are published under `<output_dir>/qc-ProteomicDDA-1.0.0/<dataset_id>/`:

| Output | Description |
| --- | --- |
| `output/<dataset_id>_QC_metrics.tsv` | Combined long-form QC table. |
| `output/mzML/<sample_name>_<fraction_id>_mzML_metrics.tsv` | Per-fraction mzML metrics. |
| `output/mzML/<sample_name>_<fraction_id>_MS1_matrix.tsv` and `.json` | Selected MS1 peaks and selection metadata. |
| `output/mzML/<sample_name>_<fraction_id>_MS1_landscape.html` | Interactive MS1 figure. |
| `output/search_idXML/<sample_name>_<search_engine>_search_idXML_metrics.tsv` | Unfiltered search metrics. |
| `output/filtered_idXML/<sample_name>_filtered_idXML_metrics.tsv` | Final filtered identification metrics. |
| `output/ProteinGroups/<dataset_id>_ProteinGroups_metrics.tsv` | Cohort inference metrics. |
| `output/ProteinQuantifier/<dataset_id>_protein_quant_metrics.tsv` | Protein abundance metrics. |
| `output/ProteinQuantifier/<dataset_id>_peptide_quant_metrics.tsv` | Peptide abundance metrics. |
| `log-*/nextflow-log/` | Nextflow trace, timeline, report, and effective parameters. |
| `log-*/process-log/` | Individual process logs. |

Each metric row records `dataset_id`, `acquisition_mode`, `layer`, `sample_name`, `fraction_id`, `search_engine`, `inference_route`, `metric_name`, `group`, `value`, `unit`, `status`, `source_file`, and `calculation_version`. Cohort protein-group rows have no sample key. Programmed LC gradient duration and retention-time error are not derived from these inputs. Longitudinal CVs and correlations require explicit repeated-QC metadata.

## Testing

The UCLA nf-test suite has LFQ ProteinProphet and TMT ProteinProphet integration cases. Set `TEST_DATA_DIR` in a local `.env` file using [`.env-template`](.env-template), then run a case inside a Slurm job:

```bash
bash test/run-nftest.sh test-integration-LFQ-ProteinProphet
```

See [test/README.md](test/README.md) for the TMT case and required test-data layout. These tests require the lab's external metapipeline outputs; they are not self-contained fixtures.

## License

Licensed under GNU General Public License version 2. See [LICENSE](LICENSE).
