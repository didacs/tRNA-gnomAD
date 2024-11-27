#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { PARSE_TRNA } from './modules/parse_trna'
include { GET_GNOMAD_DATA } from './modules/get_gnomad_data'
include { SPLIT_VCF } from './modules/split_vcf'
include { GET_CONSENSUS_FROM_VCF } from './modules/get_consensus_from_vcf'
include { MERGE_DATA } from './modules/merge_data'
include { CONCATENATE } from './modules/concatenate'

// Define parameters
params.input = "$projectDir/input.ss"
params.outdir = "results"

// Log parameters
log.info """\
    P A R S E   t R N A   P I P E L I N E
    ===================================
    input       : ${params.input}
    outdir      : ${params.outdir}
    """
    .stripIndent()

// Define the workflow
workflow {
    // Create input channel
    input_ch = channel.fromPath(params.input, checkIfExists: true)

    // Create a channel for the reference directory
    ref_dir_ch = channel.fromPath(params.ref_dir, checkIfExists: true, type: 'dir')

    // Run PARSE_TRNA process
    PARSE_TRNA(input_ch)

    // Create a new channel from the CSV output
    regions_ch = PARSE_TRNA.out.csv_output
        .splitCsv(header:true)
        .map { row -> 
            [row.id, row.chromosome, row.anticodon_start.toInteger(), row.anticodon_end.toInteger()]
        }

    // Run GET_GNOMAD_DATA process
    GET_GNOMAD_DATA(regions_ch, 
                    params.GNOMAD_release_version, 
                    params.GNOMAD_data_type, 
                    params.GNOMAD_genome_build)

    vcf_with_records_ch = GET_GNOMAD_DATA.out.vcf
        .filter { it[5] == "true" }
        .map { id, chromosome, start, end, vcf_file, has_records ->
            // Optionally remove the has_records field if not needed downstream
            [id, chromosome, start, end, vcf_file]
        }

    // Split VCF by record
    SPLIT_VCF(vcf_with_records_ch)

    split_vcf_ch = SPLIT_VCF.out.split_vcf
        .map { id, chrom, start, end, split_vcfs ->
            // Ensure split_vcfs is a list
            def processed_splits = split_vcfs instanceof List ? split_vcfs : [split_vcfs]
            [id, chrom, start, end, processed_splits]
        }
        .flatMap { id, chrom, start, end, split_vcfs ->
            split_vcfs.collect { split_vcf ->
                [id, chrom, start, end, split_vcf]
            }
        }

    // get consensus sequence
    GET_CONSENSUS_FROM_VCF(split_vcf_ch.combine(ref_dir_ch))

    MERGE_DATA(GET_CONSENSUS_FROM_VCF.out.consensus.combine(PARSE_TRNA.out.csv_output))
    // MERGE_DATA.out.toList().unique().view()
    CONCATENATE(MERGE_DATA.out.toList().unique())
}

// This block is executed when the pipeline script completes
workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}