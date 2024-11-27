process GET_CONSENSUS_FROM_VCF {
    debug true
    tag "${chrom}:${start}-${end}"
    publishDir "${params.outdir}/consensus", mode: 'copy'

    input:
    tuple val(id), val(chrom), val(start), val(end), path(vcfFile), path(ref_dir)

    output:
    tuple val(id), path(vcfFile), path("${vcfFile.simpleName}.consensus.fa"), emit: consensus
    path "versions.yml", emit: versions

    script:
    """
    # bcftools view --output-type z --output "${vcfFile.baseName}.sorted.vcf.gz" ${vcfFile}
    bcftools index "${vcfFile}"

    samtools faidx "${ref_dir}/${params.fasta_filename}" "${chrom}:${start}-${end}" \
    | bcftools consensus -f - "${vcfFile}" > "${vcfFile.simpleName}.consensus.fa"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        samtools: \$(samtools --version | head -n1 | sed 's/^.*samtools //; s/ .*\$//')
    END_VERSIONS
    """
}