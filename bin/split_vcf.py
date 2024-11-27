#!/usr/bin/env python

from pysam import VariantFile
import argparse
import os

parser = argparse.ArgumentParser(description='')
parser.add_argument('--vcf', help='', required=True)
args = parser.parse_args()

def split_vcf_by_variant(vcf):
    vcf_in = VariantFile(vcf)
    output_written = False

    for record in vcf_in:
        # Create a filename based on CHROM, POS, REF, and ALT
        chrom = record.chrom
        pos = record.pos
        ref = record.ref
        assert len(record.alts) == 1
        alt = record.alts[0]
        output_filename = f"split_{chrom}-{pos}-{ref}-{alt}.vcf.bgz"

        # Open a new VCF file for writing the single variant
        with VariantFile(output_filename, 'wb', header=vcf_in.header) as vcf_out:
            vcf_out.write(record)
        output_written = True


    if not output_written:
        output_filename = f"split_{os.path.basename(args.vcf).split('.')[0]}_empty.vcf.bgz"
        with VariantFile(output_filename, 'wb', header=vcf_in.header) as vcf_out:
            pass  # Just write the header and no records

split_vcf_by_variant(args.vcf)
