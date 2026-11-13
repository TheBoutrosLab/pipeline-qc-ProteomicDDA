process summarize_ProteinQuantifier {
    container params.docker_image_pyopenms
    containerOptions params.project_script_mount_options

    publishDir "${params.output_dir_base}/output/ProteinQuantifier", mode: 'copy'

    label 'process_low'

    input:
        tuple(
            val(level),
            path(input_file),
            path(mzml_csv, name: 'samples.csv'),
            path(channel_csv, name: 'channels.csv'),
            path(assay_mapping_csv, name: 'assays.csv')
        )

    output:
        path("${params.dataset_id}_${level}_quant_metrics.tsv"), emit: metrics

    script:
    """
    python3 ${projectDir}/src/SummarizeQuantification.py \\
        --level ${level} \\
        --input ${input_file} \\
        --dataset-id ${params.dataset_id} \\
        --quantification-mode ${params.quantification_mode} \\
        --mzml-csv samples.csv \\
        --channel-csv channels.csv \\
        --assay-mapping-csv ${params.assay_mapping_csv ? 'assays.csv' : 'none'} \\
        --output ${params.dataset_id}_${level}_quant_metrics.tsv
    """
}
