#!/usr/bin/bash
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#
# Get Ligilactobacillus murinus genome from NCBI
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

PREFIX=${1:-lm}

# Create directory
mkdir -p "$PREFIX/"

# Download genome FASTA
if [ ! -f "$PREFIX/genome.fasta" ]
then
	curl 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/288/115/GCF_003288115.1_ASM328811v1/GCF_003288115.1_ASM328811v1_genomic.fna.gz' \
	  | gzip -d > $PREFIX/genome.fasta
fi

# Download annotations GTF and GFF
# For this genome, duplicate CDS features into exon feature so mappers can understand them
if [ ! -f "$PREFIX/genome.gtf" ]
then
curl 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/288/115/GCF_003288115.1_ASM328811v1/GCF_003288115.1_ASM328811v1_genomic.gtf.gz' \
  | gzip -d \
  | awk 'BEGIN{FS=OFS="\t"}{print}($3=="CDS"){$3="exon";print}' > $PREFIX/genome.gtf
fi

if [ ! -f "$PREFIX/genome.gff.gz" ]
then
	curl 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/288/115/GCF_003288115.1_ASM328811v1/GCF_003288115.1_ASM328811v1_genomic.gff.gz' \
	  > "$PREFIX/genome.gff.gz"
fi


#-#-#-#-#-#-#-#-#-#-#
# Genome Indexing
#-#-#-#-#-#-#-#-#-#-#
[ ! -f "$PREFIX/genome.fasta.fai" ] && samtools faidx "$PREFIX/genome.fasta"

# HISAT2 genome indexing
if [ ! -d "$PREFIX/ht2_index/" ]
then
	mkdir -p "$PREFIX/ht2_index/" && hisat2-build -p 6 "$PREFIX/genome.fasta" "$PREFIX/ht2_index/index"
fi

# BWA genome indexing
if [ ! -d "$PREFIX/bwa_index/" ]
then
	mkdir -p "$PREFIX/bwa_index/" \
	  && ln -s ../genome.fasta "$PREFIX/bwa_index/" \
	  && bwa index "$PREFIX/bwa_index/genome.fasta"
fi



