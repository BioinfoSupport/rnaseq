#!/usr/bin/bash
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#
# Merged Mycobacterium marinum and Dicty genomes into one
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

PREFIX=${1:-DdMm}

# Create directory
mkdir -p "$PREFIX/"

build_ref_Dd.sh "$PREFIX/Dd"
build_ref_Mm.sh "$PREFIX/Mm"

if [ ! -f "$PREFIX/genome.fasta" ]
then
	cat "$PREFIX/Dd/genome.fasta" "$PREFIX/Mm/genome.fasta" > "$PREFIX/genome.fasta"
fi
if [ ! -f "$PREFIX/genome.gtf" ]
then
	cat "$PREFIX/Dd/genome.gtf" "$PREFIX/Mm/genome.gtf" > "$PREFIX/genome.gtf"
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



