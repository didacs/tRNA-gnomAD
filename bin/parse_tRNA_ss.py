#!/usr/bin/env python

import csv
import re
import logging
import sys
from dataclasses import dataclass, fields
from typing import List

@dataclass
class TRNAEntry:
    id: str
    chromosome: str
    start: int
    end: int
    strand: str
    length: int
    type: str
    anticodon: str
    anticodon_start: int
    anticodon_end: int
    score: float
    sequence: str
    structure: str

VALID_CHROMOSOMES = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

def parse_trnascan_output(input_file: str) -> List[TRNAEntry]:
    entries = []
    
    # Read the entire file content
    with open(input_file, 'r') as f:
        content = f.read()
    
    # Split the content into individual tRNA entries
    trna_entries = content.strip().split('\n\n')
    
    for entry in trna_entries:
        if not entry.strip():
            continue
            
        # Parse the basic information
        id_match = re.search(r'([\w.]+)\s+\((\d+)-(\d+)\)\s+Length:\s+(\d+)', entry)
        if not id_match:
            continue
            
        trna_id, start, end, length = id_match.groups()
        chromosome = trna_id.split('.')[0]

        # skip record if not in main chromosome
        if chromosome not in VALID_CHROMOSOMES:
            print(f"Skipped record on invalid chromosome: {chromosome}")
            continue
        
        # Parse type and anticodon information
        type_match = re.search(r'Type:\s+(\w+)\s+Anticodon:\s+(\w+)\s+at\s+\d+-\d+\s+\((\d+)-(\d+)\)\s+Score:\s+([\d.]+)', entry)
        if not type_match:
            continue
            
        trna_type, anticodon, ac_start, ac_end, score = type_match.groups()

        if end < start:
            start, end = end, start
            ac_start, ac_end = ac_end, ac_start
            strand = '-'
        else:
            strand = '+'
        
        # Parse sequence and structure
        seq_match = re.search(r'Seq:\s+([A-Za-z]+)', entry)
        str_match = re.search(r'Str:\s+([><.]+)', entry)
        
        sequence = seq_match.group(1) if seq_match else ''
        structure = str_match.group(1) if str_match else ''
        
        # Create TRNAEntry object
        trna_entry = TRNAEntry(
            id=trna_id,
            chromosome=chromosome,
            start=int(start),
            end=int(end),
            strand=strand,
            length=int(length),
            type=trna_type,
            anticodon=anticodon,
            anticodon_start=int(ac_start),
            anticodon_end=int(ac_end),
            score=float(score),
            sequence=sequence,
            structure=structure
        )
        entries.append(trna_entry)
    
    return entries

def write_to_csv(entries: List[TRNAEntry], output_file: str):
    with open(output_file, 'w', newline='') as f:
        # Get field names from dataclass
        fieldnames = [field.name for field in fields(TRNAEntry)]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        
        writer.writeheader()
        for entry in entries:
            writer.writerow(entry.__dict__)

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <output_file>")
        sys.exit(1)
        
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    try:
        entries = parse_trnascan_output(input_file)
        write_to_csv(entries, output_file)
        print(f"Successfully parsed {len(entries)} tRNA entries to {output_file}")
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()