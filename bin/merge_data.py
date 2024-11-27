#!/usr/bin/env python3

import argparse
import pysam
import pandas as pd
from Bio.Seq import Seq
from Bio.Data.CodonTable import standard_dna_table
import os


def translate_anticodon_to_amino_acid(anticodon: str, table: int = 1) -> str:
    """
    Translate a tRNA anticodon to its corresponding three-letter amino acid code using Biopython.
    
    Args:
        anticodon (str): The RNA anticodon sequence (e.g., "AAG").
        table (int): The translation table to use (default: 1 for the standard genetic code).
    
    Returns:
        str: Three-letter amino acid code.
    """
    # Map from single-letter to three-letter amino acid codes
    one_to_three = {
        'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu',
        'F': 'Phe', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
        'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn',
        'P': 'Pro', 'Q': 'Gln', 'R': 'Arg', 'S': 'Ser',
        'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr',
        '*': 'Stop'
    }
    
    # Reverse complement to get the corresponding codon
    anticodon_seq = Seq(anticodon)
    codon = str(anticodon_seq.reverse_complement().transcribe())
    single_letter = Seq(codon).translate(table=table)[0]
    
    # Convert to three-letter code
    return one_to_three.get(single_letter, 'Unknown')


def read_fasta(fasta_file: str) -> str:
    """
    Read sequence from a FASTA file.
    
    Args:
        fasta_file (str): Path to the FASTA file.
    
    Returns:
        str: The sequence read from the file.
    """
    with open(fasta_file, 'r') as f:
        # Skip header line
        next(f)
        # Read sequence
        sequence = next(f).strip()
    return sequence


def main(id: str, bcf_file: str, fasta_file: str, csv_file: str) -> None:
    """
    Main function to analyze tRNA variants.

    Args:
        id (str): ID of the tRNA to analyze.
        bcf_file (str): Path to the BCF/VCF file.
        fasta_file (str): Path to the FASTA file.
        csv_file (str): Path to the CSV file with tRNA details.
    """
    # Load CSV and get row for matching ID
    df = pd.read_csv(csv_file)
    row = df[df['id'] == id].iloc[0]
    trna_type: str = row['type']
    anticodon: str = row['anticodon']
    strand: str = row['strand']
    
    # Initialize variant info with default values
    variant_info: dict[str, str] = {
        'ref': 'NA',
        'alt': 'NA',
        'AF': 'NA',
        'AC': 'NA',
        'AN': 'NA',
        'nhomalt': 'NA'
    }
    
    # Try to load BCF file and get variant info if not empty
    bcf = pysam.VariantFile(bcf_file)
    try:
        record = next(bcf)
        variant_info = {
            'ref': record.ref,
            'alt': ','.join(str(alt) for alt in record.alts),
            'AF': record.info.get('AF', ['NA'])[0],
            'AC': record.info.get('AC', ['NA'])[0],
            'AN': record.info.get('AN', 'NA'),
            'nhomalt': record.info.get('nhomalt', ['NA'])[0]
        }
    except StopIteration:
        # File is empty, use default values
        pass
    
    # Load sequence from FASTA
    sequence: str = read_fasta(fasta_file)
    
    # Reverse complement if on negative strand
    if strand == '-':
        sequence = str(Seq(sequence).reverse_complement())
    
    # Translate anticodon to amino acid
    amino_acid: str = translate_anticodon_to_amino_acid(sequence)
    
    # Build results table
    results: dict[str, str] = {
        'ID': id,
        'tRNA_type': trna_type,
        'Anticodon': anticodon,
        'Strand': strand,
        'Sequence': sequence,
        'Amino_acid': amino_acid,
        'Reference': variant_info['ref'],
        'Alternative': variant_info['alt'],
        'Allele_Frequency': variant_info['AF'],
        'Allele_Count': variant_info['AC'],
        'Allele_Number': variant_info['AN'],
        'Homozygous_Alt_Count': variant_info['nhomalt']
    }
    
    # Convert to DataFrame for nice display
    results_df: pd.DataFrame = pd.DataFrame([results])
    
    # Save to file
    outfile: str = os.path.basename(bcf_file).split('.')[0]
    results_df.to_csv(f"{outfile}_merged.tsv", sep='\t', index=False)


if __name__ == "__main__":
    parser: argparse.ArgumentParser = argparse.ArgumentParser(description='Analyze tRNA variants')
    parser.add_argument('--id', required=True, help='ID to analyze')
    parser.add_argument('--bcf', required=True, help='BCF/VCF file path')
    parser.add_argument('--fasta', required=True, help='FASTA file path')
    parser.add_argument('--csv', required=True, help='CSV file path')
    
    args: argparse.Namespace = parser.parse_args()
    main(args.id, args.bcf, args.fasta, args.csv)
