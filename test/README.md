# Integration tests

Run either ProteinProphet case inside a Slurm allocation:

- `bash test/run-nftest.sh test-integration-LFQ-ProteinProphet`
- `bash test/run-nftest.sh test-integration-TMT-ProteinProphet`

Both cases use repository `mzml.csv`, `search.csv`, `filtered.csv`, and `features.csv` manifests with explicit test-data paths. TMT also uses `channel.csv` with the biological reporter-channel descriptions. The feature manifests point to published PSMFeatureExtractor idXML files from validation 2.3.2. Copy `.env-template` to `.env` and set `TEST_DATA_DIR`. The CSV path references are expanded from that variable at test launch.

Each UCLA nf-test case checks the merged QC table for representative acquisition, identification, protein-inference, and quantification metrics, including TMT biological `group` and reporter-channel identity.

The LFQ case uses six mzML fractions across two samples from `$TEST_DATA_DIR/ecoli_lfq/metapipeline-ProteomicDDA/metapipeline-ProteomicDDA-1.2.7`. The TMT case uses six mzML fractions across two plexes from `$TEST_DATA_DIR/ecoli_tmt/metapipeline-ProteomicDDA_FeatureBased/metapipeline-ProteomicDDA-1.2.7`. Both use current ProteinQuantifier 1.3.0 output with `abundance|...|ch*` headers and the published filtered ProteinProphet group idXML.
