# Integration tests

Run either ProteinProphet case inside a Slurm allocation:

- `bash test/run-nftest.sh test-integration-LFQ-ProteinProphet`
- `bash test/run-nftest.sh test-integration-TMT-ProteinProphet`

Both cases use repository `mzml.csv`, `search.csv`, `filtered.csv`, and `features.csv` manifests with explicit test-data paths. The feature manifests point to published PSMFeatureExtractor idXML files from validation 2.3.1. Copy `.env-template` to `.env` and set `TEST_DATA_DIR`. The CSV path references are expanded from that variable at test launch.

The LFQ case uses four mzML fractions across two samples from `$TEST_DATA_DIR/ecoli_lfq/metapipeline-ProteomicDDA-ProteinProphet`. The TMT case uses six mzML fractions across two plexes from the corresponding `ecoli_tmt` run. Both use metapipeline 1.0.1 outputs, including the published filtered ProteinProphet group idXML from ProteinQuantifier 1.1.2.
