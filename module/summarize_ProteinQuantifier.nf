process summarize_ProteinQuantifier {
    container params.docker_image_pyopenms
    containerOptions params.project_script_mount_options
    publishDir "${params.output_dir_base}/output/ProteinQuantifier", mode: 'copy'
    label 'process_low'
    input:
        tuple val(level), path(input_file)
    output:
        path("${params.dataset_id}_${level}_quant_metrics.tsv"), emit: metrics
    script:
    """
    python3 ${projectDir}/src/SummarizeQuantification.py --level ${level} --input ${input_file} --dataset-id ${params.dataset_id} --output ${params.dataset_id}_${level}_quant_metrics.tsv
    """
}
