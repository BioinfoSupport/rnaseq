#!/usr/bin/bash
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#
# Get Dicty genome from NCBI and concatenate them
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

PREFIX=${1:-Dd}

# Create directory
mkdir -p "$PREFIX/"

# Download genome FASTA
if [ ! -f "$PREFIX/genome.fasta" ]
then
	curl 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/695/GCF_000004695.1_dicty_2.7/GCF_000004695.1_dicty_2.7_genomic.fna.gz' \
	  | gzip -d > "$PREFIX/genome.fasta"
fi

# Download annotations GTF
if [ ! -f "$PREFIX/genome.gtf" ]
then
	curl 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/695/GCF_000004695.1_dicty_2.7/GCF_000004695.1_dicty_2.7_genomic.gtf.gz' \
	  | gzip -d > "$PREFIX/genome.gtf"
fi

if [ ! -f "$PREFIX/genome.gff.gz" ]
then
	curl 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/695/GCF_000004695.1_dicty_2.7/GCF_000004695.1_dicty_2.7_genomic.gff.gz' \
	  > "$PREFIX/genome.gff.gz"
fi

#-#-#-#-#-#-#-#-#-#-#
# Genome Indexing
#-#-#-#-#-#-#-#-#-#-#
[ ! -f "$PREFIX/genome.fasta.fai" ] && samtools faidx "$PREFIX/genome.fasta"

# HISAT2 genome indexing
if [ ! -d "$PREFIX/ht2_index/" ]
then
	mkdir -p "$PREFIX/ht2_index/" \
	  && hisat2_extract_exons.py "$PREFIX/genome.gtf" > "$PREFIX/ht2_index/genome.exon" \
	  && hisat2_extract_splice_sites.py "$PREFIX/genome.gtf" > "$PREFIX/ht2_index/genome.ss" \
	  && hisat2-build -p 6 --exon "$PREFIX/ht2_index/genome.exon" --ss "$PREFIX/ht2_index/genome.ss" "$PREFIX/genome.fasta" "$PREFIX/ht2_index/index"
fi

# BWA genome indexing
if [ ! -d "$PREFIX/bwa_index/" ]
then
	mkdir -p "$PREFIX/bwa_index/" \
	  && ln -s ../genome.fasta "$PREFIX/bwa_index/" \
	  && bwa index "$PREFIX/bwa_index/genome.fasta"
fi



