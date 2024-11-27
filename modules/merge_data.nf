process MERGE_DATA {
    publishDir "results/merge_data", mode: 'copy'
    
    input:
        tuple val(id), path(vcf_file), path(fasta), path(csv)
    
    output:
        path "${vcf_file.simpleName}_merged.tsv"
    
    script:
    """
    bcftools view --output-type v --output ${vcf_file}.uncompressed.vcf ${vcf_file};
    merge_data.py \
        --id "${id}" \
        --bcf ${vcf_file}.uncompressed.vcf \
        --fasta "${fasta}" \
        --csv "${csv}"
    """
}