process CONCATENATE {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(tsv_files)

    output:
    path "merged_table.tsv", emit: merged_table

    script:
    """
    #!/usr/bin/env bash
    set -e
    set -x

    # Check if there are any files
    echo "Number of input files: ${tsv_files.size()}"
    if [ ${tsv_files.size()} -eq 0 ]; then
        echo "No files to merge" >&2
        exit 1
    fi

    # If only one file, just copy it
    if [ ${tsv_files.size()} -eq 1 ]; then
        cp ${tsv_files[0]} merged_table.tsv
        echo "Copied single file to merged_table.tsv"
        exit 0
    fi

    # Multiple files: merge with header from first file
    first_file=${tsv_files[0]}
    
    # Print header
    head -n 1 "\$first_file" > merged_table.tsv
    echo "Added header from \$first_file to merged_table.tsv"
    
    # Append data rows from all files, skipping their headers
    for file in ${tsv_files}; do
        echo "Processing file: \$file"
        tail -n +2 "\$file" >> merged_table.tsv
    done

    # Check the size of the output file
    echo "Size of merged_table.tsv: \$(wc -l < merged_table.tsv) lines"

    # Display the first few lines of the output file
    echo "First few lines of merged_table.tsv:"
    head -n 5 merged_table.tsv
    """
}