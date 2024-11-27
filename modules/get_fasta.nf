process GET_FASTA {
    tag "$bed"
    publishDir "${params.outdir}/fasta", mode: 'copy'

    input:
    path bed
    path ref_dir

    output:
    path "${bed.baseName}.fasta", emit: fasta
    path "versions.yml", emit: versions

    script:
    """
    bedtools getfasta -s -fi ${ref_dir}/${params.fasta_filename} -bed $bed -fo ${bed.baseName}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}