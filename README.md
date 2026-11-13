# pipeline-qc-ProteomicDDA

This Nextflow DSL2 pipeline collects QC metrics from completed `metapipeline-ProteomicDDA` output. It reads each original mzML fraction, unfiltered search idXML for each sample and engine, the final filtered validation idXML for each sample or TMT plex, optional published PSM feature files, one cohort protein-inference artifact, and the ProteinQuantifier peptide and protein tables. It also creates an interactive 3D MS1 peak map for every mzML fraction.

![DDA QC workflow schematic](img/flow-diagram.svg?raw=true)

## How to run

1. Initialize Git submodules with `git submodule update --init --recursive`.
2. Copy [config/template.config](config/template.config) to a run-specific config and set all input paths, `quantification_mode`, `inference_route`, `dataset_id`, and `output_dir`.
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

The search CSV can be the metapipeline's published `validation.csv`, and the filtered CSV can be its published `quant-idxml.csv`. For TMT, `sample_name` is the plex ID in the mzML and search manifests, while the filtered manifest uses `plex_id,idXML` as published by the metapipeline. LFQ uses `sample_name,idXML` in the filtered manifest. Each mzML sample or plex must occur in both identification manifests. Use one final consensus filtered idXML per sample when searches from multiple engines were merged; do not sum it with per-engine filtered counts.

Set `quantification_mode` explicitly to `TMT` or `LFQ`. TMT also requires `channel_csv` with `plex_id,channel,description`, matching the reporter-channel manifest passed to TMT quantitation. The `description` is the biological sample or reference label written to quantification metric `group`. Unassigned channels have an empty `group`.

If a ProteinQuantifier abundance header cannot be matched unambiguously to an mzML sample/plex, provide `assay_mapping_csv` with columns `level,abundance_column,sample_name,fraction_id,reporter_channel,group`. `level` is `peptide` or `protein`; `abundance_column` is the exact table header. Leave `fraction_id` blank for an assay spanning fractions, and leave `group` blank for an unassigned TMT channel. The pipeline fails on ambiguous headers rather than guessing from a filename or fraction-group number.

**Optional PSM feature CSV** (`psm_features_csv`):

```csv
sample_name,search_engine,route,fraction_id,path
sample_01,MSFragger,MSBooster,1,/data/sample_01/MSBooster/output/engine_MSFragger/run_01-fraction_1_edited.pin
sample_01,Comet,PSMFeatureExtractor,,/data/sample_01/OpenMS-3.5.0/output/sample_01_Comet_featureExtracted.idXML
```

`route` identifies the producing process: `MSBooster` for an edited PIN or `PSMFeatureExtractor` for an idXML. List one row per published file. MSBooster rows require a `fraction_id` present in the mzML CSV; PSMFeatureExtractor rows may leave it blank. The sample and search engine must occur together in the search CSV. A header-only CSV is valid when no feature files were published. Paths may be absolute or relative to the feature CSV.

| Parameter | Description |
| --- | --- |
| `dataset_id` | Cohort identifier used in the output path and metric filenames. |
| `quantification_mode` | Required `TMT` or `LFQ` mode. |
| `channel_csv` | TMT plex/channel/description CSV; required for TMT. |
| `assay_mapping_csv` | Optional explicit mapping for ambiguous abundance headers. |
| `mzml_csv` | Per-fraction mzML manifest. |
| `search_idxml_csv` | Unfiltered search idXML manifest, keyed by sample and search engine. |
| `filtered_idxml_csv` | Final filtered idXML manifest, keyed by sample or TMT plex. |
| `psm_features_csv` | Optional manifest of published MSBooster PIN or PSMFeatureExtractor idXML files. |
| `inference_route` | `ProteinProphet` or `ProteinInferencer`. |
| `protein_groups` | Cohort filtered ProteinProphet idXML or ProteinInferencer consensusXML. |
| `protein_quant_csv` | ProteinQuantifier protein abundance table. |
| `peptide_quant_csv` | ProteinQuantifier peptide abundance table. |
| `output_dir` | Root directory for QC results. |
| `work_dir` | Optional Nextflow work directory. |
| `pyopenms_version` | pyOpenMS image tag; default `3.5.0`. |
| `docker_image_plotly` | Plotly image; default `multiqc/multiqc:v1.35`. |
| `base_resource_update` | Optional process CPU and memory overrides. |

ProteinProphet grouping uses the filtered idXML passed to ProteinQuantifier, now published in the protein quantitation output directory. ProteinInferencer grouping uses the inferred consensusXML; choose the protein-FDR-filtered version if that was passed downstream. The two current ProteinQuantifier tables are tab-delimited despite their `.csv` extension and use `abundance|...|ch*` headers. QC uses those headers to map each abundance column to a sample or TMT plex; TMT quantification uses the biological channel description in `group`, while LFQ uses `group=all`. Older ProteinQuantifier header layouts are not accepted.

## Pipeline steps

| Process | Input | Output |
| --- | --- | --- |
| `summarize_mzML` | One original mzML fraction | Acquisition duration, MS1/MS2 scan counts and rates, survey cycle timing, MS2 scans per cycle, and whole-run TIC summary |
| `extract_MS1_matrix` | One mzML fraction | Bounded RT/m/z/intensity peak matrix and selection metadata |
| `plot_MS1_landscape` | Matrix and metadata | Offline interactive 3D MS1 stick spectrum |
| `summarize_search_idXML` | One raw search idXML per sample and engine | Candidate and top-hit counts, distinct top-hit peptides, score metadata, unknown target/decoy labels, signed/absolute ppm mass error, and precursor charge summaries |
| `summarize_filtered_idXML` | One final filtered idXML per sample or plex | All accepted PSMs, distinct retained peptides, validated score metadata, and unknown target/decoy labels |
| `summarize_PSMFeatures` | One published MSBooster PIN or PSMFeatureExtractor idXML | Target/decoy feature counts and distributions, including mass and RT deviations when available |
| `summarize_ProteinGroups` | One cohort inference artifact | Group count and size, plus ProteinProphet probability or ProteinInferencer score summaries when present |
| `summarize_ProteinQuantifier` | Peptide and protein-group abundance tables | Distinct quantified peptide/protein-group counts and log10 positive intensity summaries per assay column |
| `collect_QCMetrics` | All stage TSVs | Combined `<dataset_id>_QC_metrics.tsv` |

A continuous distribution is saved as its nonmissing count, missing count, and 5th, 25th, 50th, 75th, and 95th percentiles. `group=all` includes target, decoy, and unknown labels; these classes also have separate summaries. A required distribution whose source field is absent has all seven rows marked `unavailable`. Raw top hits mean the best finite-score hit per idXML peptide-identification record, not necessarily one hit per spectrum. Charge states have counts and fractions. A recorded mzML acquisition duration is not necessarily the programmed LC gradient duration.

PSM feature summaries use actual PIN columns or OpenMS `extra_features` metadata: `deltamass`, `dm`, `absdm`, RT delta or difference fields, and spectral angle, entropy, correlation, or cosine fields when present. `deltamass` is an m/z difference; when `CalcMass` is present, signed and absolute ppm mass errors are also reported. Other source feature units are left blank unless established by the source. Missing derived ppm mass-error families are marked `unavailable` rather than zero; absent optional RT and spectral feature columns emit no rows. Counts and distributions are reported for all PSMs and separately for targets, decoys, and unknown labels.

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
| `output/PSMFeatures/<sample_name>_<search_engine>_<route>_<fraction_id>_PSMFeatures_metrics.tsv` | Published feature metrics, when provided. |
| `output/ProteinGroups/<dataset_id>_ProteinGroups_metrics.tsv` | Cohort inference metrics. |
| `output/ProteinQuantifier/<dataset_id>_protein_quant_metrics.tsv` | Protein abundance metrics. |
| `output/ProteinQuantifier/<dataset_id>_peptide_quant_metrics.tsv` | Peptide abundance metrics. |
| `log-*/nextflow-log/` | Nextflow trace, timeline, report, and effective parameters. |
| `log-*/process-log/` | Individual process logs. |

Each metric row records the same 12 columns, in the same order, as the DIA QC pipeline: `dataset_id`, `layer`, `source_process`, `sample_name`, `search_engine`, `fraction_id`, `metric_name`, `group`, `value`, `unit`, `status`, `source_file`. The merged table follows this layer order: `acquisition`, `raw_identification`, `psm_features`, `validated_identification`, `protein_inference`, `spectral_library`, `quantification`. The current DDA pipeline has no spectral-library metrics, so that layer has no rows. Cohort protein-group rows have no sample key. Programmed LC gradient duration is not derived from these inputs. RT deviation requires a published prediction feature. Longitudinal CVs and correlations require explicit repeated-QC metadata. This schema is a breaking change: no old-name aliases or legacy parser are provided.

## Testing

The UCLA nf-test suite has LFQ ProteinProphet and TMT ProteinProphet integration cases. Set `TEST_DATA_DIR` in a local `.env` file using [`.env-template`](.env-template), then run a case inside a Slurm job:

```bash
bash test/run-nftest.sh test-integration-LFQ-ProteinProphet
```

See [test/README.md](test/README.md) for the TMT case and required test-data layout. These tests require the lab's external metapipeline outputs; they are not self-contained fixtures.

## License

Licensed under GNU General Public License version 2. See [LICENSE](LICENSE).
