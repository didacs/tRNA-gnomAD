# tRNA anticodon genetic variants
This pipeline queries gnomAD SNVs that impact the anticodon in tRNAs. Such variants could potentially have an effect on mRNA translation as they create a mismatch between the codon the tRNA recognizes and the amino acid it inserts, essentially misreading a correct codon in the CDS.

## install nextflow
Follow steps in https://www.nextflow.io/docs/latest/install.html#install-nextflow

## Create conda env
Not necessary if samtools is not in your PATH, which is required to index the genome fasta file below
```
mamba env create -f environment.yml
mamba activate trna-vars
```

## Download and index hg38 reference genome
```
mkdir ref
curl -o ref/GRCh38.p14.genome.fa.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.p14.genome.fa.gz
gunzip ref/GRCh38.p14.genome.fa.gz
samtools faidx ref/GRCh38.p14.genome.fa
```

## Download genomic tRNA set from gtRNAdb
```
curl -O https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.tar.gz
mkdir hg38-tRNAs
tar -xvzf hg38-tRNAs.tar.gz -C hg38-tRNAs
```

## run the pipeline
```
nextflow run /Users/didac/github/tRNA-gnomAD/main.nf --input hg38-tRNAs/hg38-tRNAs-confidence-set.ss --ref_dir ref/ --fasta_filename GRCh38.p14.genome.fa
```

