# Integration tests

Run either ProteinProphet case inside a Slurm allocation:

- `bash test/run-nftest.sh test-integration-LFQ-ProteinProphet`
- `bash test/run-nftest.sh test-integration-TMT-ProteinProphet`

Both cases use repository `mzml.csv`, `search.csv`, and `filtered.csv` manifests with explicit test-data paths. Copy `.env-template` to `.env` and set `TEST_DATA_DIR`. The CSV path references are expanded from that variable at test launch.

The LFQ case uses four mzML fractions across two samples from `$TEST_DATA_DIR/ecoli_lfq/metapipeline-ProteomicDDA-ProteinProphet`. The TMT case uses six mzML fractions across two plexes from the corresponding `ecoli_tmt` run. Both cases require the filtered ProteinProphet group idXML used by ProteinQuantifier, which is published when the upstream protein quantitation run uses `save_intermediate_files = true`.
