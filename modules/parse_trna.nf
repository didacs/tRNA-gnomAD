process PARSE_TRNA {
    tag "$input_file.simpleName"
    publishDir "${params.outdir}/parsed_trna", mode: 'copy'

    input:
    path input_file

    output:
    path "${input_file.simpleName}.csv", emit: csv_output
    path "versions.yml", emit: versions

    script:
    """
    parse_tRNA_ss.py $input_file ${input_file.simpleName}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}