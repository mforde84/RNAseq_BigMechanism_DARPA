# DARPA_UChicago_RNAseq_pipeline

usage ./pipeline.sh genome_directory fastq_directory

output is in $(pwd)/output

Dependancies are based upon Ubuntu LTS 14.04. 
Different environments will likely require scripting modifications.

Genome directory requires a merged fasta primary assembly (e.g., ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz) and the associated GTF annotation file (e.g., ftp://ftp.ensembl.org/pub/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh38.87.gtf.gz)
