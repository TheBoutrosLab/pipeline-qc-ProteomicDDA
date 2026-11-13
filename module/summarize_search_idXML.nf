process summarize_search_idXML {
    container params.docker_image_pyopenms
    containerOptions params.project_script_mount_options

    publishDir "${params.output_dir_base}/output/search_idXML", mode: 'copy'

    label 'process_low'

    input:
        tuple(
            val(sample_name),
            val(search_engine),
            path(input_file)
        )

    output:
        path("${sample_name}_${search_engine}_search_idXML_metrics.tsv"), emit: metrics

    script:
    """
    python3 ${projectDir}/src/SummarizeIdXML.py \\
        --stage search \\
        --sample-name ${sample_name} \\
        --search-engine ${search_engine} \\
        --input ${input_file} \\
        --dataset-id ${params.dataset_id} \\
        --output ${sample_name}_${search_engine}_search_idXML_metrics.tsv
    """
}
