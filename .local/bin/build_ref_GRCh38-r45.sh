#!/usr/bin/bash
# Build human genome from gencode

PREFIX=${1:-GRCh38-r45}

# Create directory
mkdir -p "$PREFIX/"


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Download genome FASTA and annotations
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# FASTA
if [ ! -f "$PREFIX/genome.fasta" ]
then
	curl 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh38.p14.genome.fa.gz' \
	  | gzip -d > "$PREFIX/genome.fasta"
fi

# GTF
if [ ! -f "$PREFIX/genome.gtf" ]
then
	curl 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz' \
	  | gzip -d > "$PREFIX/genome.gtf"
fi

# GFF3
if [ ! -f "$PREFIX/genome.gff3.gz" ]
then
	curl 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gff3.gz' \
	  > "$PREFIX/genome.gff3.gz"
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

