#!/usr/bin/bash
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#
# Merged Human and lm genomes
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

PREFIX=${1:-GRCh38-lm}

# Create directory
mkdir -p "$PREFIX/"

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Make each genome individually
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
build_ref_GRCh38-r45.sh "$PREFIX/GRCh38-r45"
build_ref_lm.sh "$PREFIX/lm"

#-#-#-#-#-#-#-#-#-#-#
# Merge them
#-#-#-#-#-#-#-#-#-#-#
if [ ! -f "$PREFIX/genome.fasta" ]
then
	cat "$PREFIX/GRCh38-r45/genome.fasta" "$PREFIX/lm/genome.fasta" > "$PREFIX/genome.fasta"
fi
if [ ! -f "$PREFIX/genome.gtf" ]
then
	cat "$PREFIX/GRCh38-r45/genome.gtf" "$PREFIX/lm/genome.gtf" > "$PREFIX/genome.gtf"
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



