



# 1- Initialize a new project
git clone 'https://github.com/BioinfoSupport/rnaseq-pipeline.git' my_new_project
cd my_new_project

# 2- (Build) and Run the container 
docker compose up -d

# 3- Connect to the container with a web-browser
open 'http://localhost:8787'

# 4- Run QC, Hisat2 mapping, and feature counting
rnaseq data/fastq/tests/ALL

# 5- Aggregate feature counts


# Stop the container
docker compose down



# Simplified procedure:
docker compose up -d

# Map and quantify all .fastq.gz files located in data/fastq/pilot
docker compose exec rnaseq make test
docker compose exec rnaseq make data/fastq/pilot/ALL

# Get a Bash in the container
docker compose exec rnaseq bash




docker compose exec rnaseq build_ref_DdMm.sh data/ref/DdMm
docker compose exec rnaseq build_ref_GRCh38-r45.sh data/ref/GRCh38-r45
docker compose exec rnaseq build_ref_GRCm39-M34.sh data/ref/GRCm39-M34


docker run --rm -it -v "$PWD:/home/rstudio/workdir" unigebsp/rnaseq-pipeline bash -c 'build_ref_DdMm.sh /home/rstudio/workdir/data/ref/DdMm'


# 2- If necessary, update DdMm reference genome in the skeleton
docker run --rm -it -v "$PWD:/home/rstudio/workdir" unigebsp/rnaseq-pipeline bash -c 'build_ref_DdMm.sh /home/rstudio/workdir/data/ref/DdMm'



