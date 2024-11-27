process SPLIT_VCF {
   input:
        tuple val(id), val(chrom), val(start), val(end), path(vcfFile)

    output:
        tuple val(id), val(chrom), val(start), val(end), path("split_*.vcf.bgz"), emit: split_vcf
    script:
    """
    bcftools index ${vcfFile};
    split_vcf.py \
        --vcf ${vcfFile}
    """
}