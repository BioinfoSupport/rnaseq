


# About
This repository implements an RNA-seq pipeline, which is able to:

 1. map reads in FASTQ.GZ format with HISAT2
 2. quantify reads in a GTF file
 3. produce Quality Control plots with R

All required software are packaged into a (docker-)container which make it runnable on any platform.



# Quick start

 1) Make sure `docker` (https://www.docker.com) is installed on your system and is running.
 
 2) Download or clone the current repository.
```
git clone 'https://github.com/BioinfoSupport/rnaseq-pipeline.git' my_new_project
```

 3) Copy your `.fastq.gz` files into subfolder `data/fastq/`.
 
 4) Put your reference genome into subfolder `data/ref`.
 
 5) Run the container
```
docker compose up -d
```

  6.1) Launch the pipeline within the container
```
docker compose exec rnaseq GENOME=DdMm data/fastq/ALL
```
	6.2) Alternatively connect to the Rstudio GUI running at URL http://localhost:8787



# Prepare a reference genome

A reference genome has to be provided to the pipeline to work properly, and has to be prepared before any read alignment can take place. The reference genome folder you provide must follow a this specific structure:
```
data/ref/DdMm
  genome.fasta       # DNA sequence (downloaded from NCBI)
  genome.gtf         # Genes annotations (downloaded from NCBI)
  genome.fasta.fai   # Samtools index of the sequence (generated from FASTA)
  ht2_index/         # HISAT2 index of the genome (generated from FASTA+GTF)
```
But each genome must be made of 4 elements: the reference sequence (`genome.fasta`) and genes annotations (`genome.gtf`) are typically downloaded from NCBI (e.g. for [_S. aureus_](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000013425.1/)); 
the `.fai` file is generated with the software `samtools` (`samtools faidx genome.fasta`);
the HISAT2 indexed directory (`ht2_index`) is generated from the fasta and the gtf with HISAT2 software.

Example scripts are provided in `.local/bin/build_ref_*.sh` to show how to generate a reference genomes folder. If you are working on _Homo Sapiens_, _Mus Musculus_, _Dictyostelium discoideum_ or _Mycobacterium marinum_ you can directly run one of the script:
```
build_ref_DdMm.sh data/ref/DdMm
build_ref_GRCh38-r45.sh data/ref/GRCh38-r45
build_ref_GRCm39-M34.sh data/ref/GRCm39-M34
```


# Directory structure

```
data/
  ref/      folder containing reference genome subdirectories
    DdMm/   a reference genome folder for _Dictyostelium discoideum_ and _Mycobacterium marinum_
  fastq/    folder containing sequenced reads (.fastq.gz)
    test/   example reads
src/
  00_data.R
  01_qc.R
.local/     hidden folder with pipeline scripts
```



# Running on a HPC cluster
```
singularity exec 'docker://unigebsp/ngs' bash
```


# TODO

## Remove the container and have the MAkefile point to bin folder directly
## Show how to run on Baobab ?
## Add samples metadata 
## Make more clear how to make a new project ()




