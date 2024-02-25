


# About
This repository implements an RNA-seq pipeline, which is able to:

 1. map reads in `.fastq.gz` format with HISAT2
 2. quantify gene expression
 3. produce quality control plots with R

This pipeline should be run in our [`ngs`](https://github.com/BioinfoSupport/ngs) container 


# Quick start

 1) Make sure `docker` (https://www.docker.com) is installed on your system and is running.
 
 2) Download or clone the repository containing the pipeline.
```
git clone 'https://github.com/BioinfoSupport/rnaseq.git' my_new_project
```

 3) Copy your `.fastq.gz` files into subfolder `data/fastq/`. Note, you can download all FASTQ files sequenced by iGE3 genomic platform from the given link with:
```
wget -P ./data/fastq --content-disposition --trust-server-names -i 'https://data.ige3.genomics.unige.ch/dataset/download/xxxxxxx.txt'
```

 4) [Optional] Put your reference genome folder into `data/ref`. Example genome folders can be found in our [`genomes` repository](https://github.com/BioinfoSupport/genomes/releases). If you are using a known reference genome, it will be downloaded automatically from the repository.
 
 5) Run the container: `docker compose up -d`

 6) Launch the pipeline within the container: `docker compose exec rnaseq make GENOME=Dd+Mm data/fastq/ALL`

 7) Connect to the Rstudio GUI running at URL http://localhost:8787 and run notebooks in `src/`




# Directory structure

```
data/
  ref/      folder containing reference genome subdirectories
    Dd+Mm/   a reference genome folder for _Dictyostelium discoideum_ and _Mycobacterium marinum_
  fastq/    folder containing sequenced reads (.fastq.gz)
    test/   example reads
    
src/
  00_qc_genome.Rmd   Notebook to compute statistics on a reference genome
  01_qc_mapping.Rmd  Notebook to extract mapping statistics
  02_DESeq.Rmd       Notebook with an example usage of the pipeline with DESeq2
  
.local/     hidden folder with pipeline-specific scripts
```


# Useful commands
```
# Run the container 
docker compose up -d

# Run tests (requires reference genomes Dd,Mm,Dd+Mm)
docker compose exec rnaseq make test

# Map and quantify all .fastq.gz files located in data/fastq/pilot
docker compose exec rnaseq make data/fastq/pilot/ALL

# Get a Bash in the container
docker compose exec rnaseq bash

# Stop the container
docker compose down
```


# Running the pipeline on a HPC cluster
```
singularity exec 'docker://unigebsp/ngs' make GENOME=Dd+Mm data/fastq/ALL
```






