
This `data/ref/` directory should contains your reference genomes and must 
follow this specific structure:
```
data/ref/
  DdMm.fasta             # DNA sequence (downloaded from NCBI)
  DdMm.gtf               # Genes annotations (downloaded from NCBI)
  DdMm.fasta.fai         # Samtools index of the sequence (generated from FASTA)
  DdMm.ht2_index/       # STAR index of the genome (generated from FASTA+GTF)
```
The folder can contain multiple reference genomes with distinct identifiers. But each genome must be made of 4 elements: the reference sequence (extension `.fasta`) and genes annotations (extension `.gtf`) are typically downloaded from NCBI (e.g. for [_S. aureus_](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000013425.1/)); 
the `.fai` file is generated with the software `samtools`;
the HISAT2 indexed directory (extension `.ht2_index`) is generated from the fasta and the gtf with STAR software. 

Example scripts are provided in `./local/bin/build_ref_*.sh` to show how to generate a reference genomes folder. If you are working on _Homo Sapiens_, _Mus Musculus_ or _Dictyostelium discoideum_ + _Mycobacterium marinum_ you can directly run one of the script:
```
mkdir -p data/ref
cd data/ref
build_ref_GRCh38-v41.sh
build_ref_GRCm39-vM27.sh
build_ref_DdMm.sh
```

Note `xxx.fasta.fai` and `xxx.star_index/` can be generated automatically if `xxx.fasta` and `xxx.gtf` exist in the folder using command:
```
rnaseq xxx.fasta.fai xxx.star_index/
```

