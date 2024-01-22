

# About
This repository implements an RNA-seq pipeline, which is able to:

 1. map reads contained into a FASTQ using STAR
 2. quantify reads in genes 
 3. produce Quality Control plots with R

All required software are packaged into a (docker-)container which make it runnable on any platform.



# Installation

Make sure `docker` (https://www.docker.com) is installed on your system and is running.



# Quick start

1. Enter the container (using the pre-build image), and mount a working directory on `/home/rstudio/workdir`. 

## {.tabset}

### Docker GUI
![Docker parameters](docker_param_screenshot.png)

### MacOS Shell
```bash
docker run --rm -p 8787:8787 -e DISABLE_AUTH=true -v "$HOME/Documents/docker/:/home/rstudio/workdir" unigebsp/rnaseq-pipeline
```
### Windows Shell
```bash
docker run --rm -p 8787:8787 -e DISABLE_AUTH=true -v "/c/Users/jahn_/Documents/docker/:/home/rstudio/workdir" unigebsp/rnaseq-pipeline
```

## {-}


2. Copy an existing reference genome to your working directory.
We recommend to put the reference genomes into the subfolder `data/ref`. 
If the reference genome is not available, see the section below to prepare a new reference genome.

3. Copy your `.fastq.gz` files containing your sequenced reads into your working directory. 
We recommend you to put them into the subfolder `data/fastq`.

4. At this step your working directory should look like this:
```bash
find /home/rstudio/workdir/
# data/ref/
#   DdMm.fasta
#   DdMm.fasta.fai
#   DdMm.gtf
#   DdMm.star_index/
#
# data/fastq/
#   batch1/
#   batch1/sample1_100k.fastq.gz
#   batch1/sample2_100k.fastq.gz
```


5. Map the `.fastq.gz` files found in `data/fastq` onto the reference genome.
```bash
rnaseq GENOME=DdMm data/fastq/STAR
```

6. Get the R notebooks and run them
```bash
cp -r /pipeline/local/scripts /workdir
```

XXX. TEST
```bash
R 
source("../repository/local/scripts/00_quantify.R")
rmarkdown::render("../repository/local/scripts/01_qc.Rmd",knit_root_dir=normalizePath("."),output_dir=normalizePath("."))

```



# Rstudio GUI
The container can also run under Rstudio GUI:
Then open your browser at URL http://localhost:8787

# Prepare a reference genome

A reference genome has to be provided to the pipeline to work properly, and has to be prepared before any read alignment can take place.
The reference genome folder you provide must follow a this specific structure:
```
data/ref/
  DdMm.fasta             # DNA sequence (downloaded from NCBI)
  DdMm.gtf               # Genes annotations (downloaded from NCBI)
  DdMm.fasta.fai         # Samtools index of the sequence (generated from FASTA)
  DdMm.star_index/       # STAR index of the genome (generated from FASTA+GTF)
```
The folder can contain multiple reference genomes with distinct identifiers. But each genome must be made of 4 elements: the reference sequence (extension `.fasta`) and genes annotations (extension `.gtf`) are typically downloaded from NCBI (e.g. for [_S. aureus_](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000013425.1/)); 
the `.fai` file is generated with the software `samtools`;
the STAR indexed directory (extension `.star_index`) is generated from the fasta and the gtf with STAR software. 

Example scripts are provided in `local/bin/build_ref_*.sh` to show how to generate a reference genomes folder. If you are working on _Homo Sapiens_, _Mus Musculus_ or _Dictyostelium discoideum_ + _Mycobacterium marinum_ you can directly run one of the script:
```
mkdir -p data/ref
cd data/ref
build_ref_GRCh38-v41.sh
build_ref_GRCm39-vM27.sh
build_ref_DdMm.sh
```

Otherwise we have implemented rules to facilitated the process: if `x.fasta` and `x.gtf` exist in the folder, you can generate `x.fasta.fai` and `x.star_index/` with the command

```bash
rnaseq x.fasta.fai x.star_index/ 
```


# Build container image (advanced users)

You can build the docker image yourself from source with the command:
```bash
docker buildx build --push --platform linux/arm64,linux/amd64 -t unigebsp/rnaseq-pipeline ./
```

Run a bash 
```bash
docker run --rm -it -v "$HOME/Documents/docker/:/home/rstudio/workdir" unigebsp/rnaseq-pipeline bash
```







# Directory structure

```
data/
  ref/    folder containing reference genome (.fasta and .gtf)
  fastq/  folder containing sequenced reads (.fastq.gz)
scripts/
  00_data.R
  01_qc.R
Makefile - file containing pipeline.
```






# TODO

## Rebuild the container with BWA+Hisat installation
## Simplify quantification step
## Try to use HiSAT2 ?
## Change location of STAR temporary files to speed up process ?
## Show how to run on Baobab ?
## Add samples metadata 
## Make more clear how to make a new project ()

