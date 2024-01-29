
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Parameters
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
NJOB = 3
NCPU = 6
GENOMEDIR = ./data/ref
GENOME = DdMm


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Set environment variables
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
MKDIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
export LOCALDIR ?= $(MKDIR)/.local/
export PATH := $(LOCALDIR)/bin/:$(PATH)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Help rule
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
usage:
	@echo "usage: rnaseq <option> target"
	@echo "options:" 
	@echo "    -k                Continue on error"
	@echo "    -j <int>          Parallel execution on <int> jobs"
	@echo "    -n                Dry run: don't run anything, just show commands"
	@echo "    GENOMEDIR='path'  Directory containing the reference genomes [$(GENOMEDIR)]"
	@echo "    GENOME='path'     Default reference genome [$(GENOME)]"
	@echo "targets:"
	@echo "    <dir>/ALL               Run <dir>/FASTQC and <dir>/HT2"
	@echo "    <dir>/FASTQC            Run fastqc on all .fastq.gz files in the directory (recursively)"
	@echo "    <dir>/HT2               Map all .fastq.gz in the directory on the default genome (recursively)"	
	@echo "    fqfile.HT2              Map a single .fastq.gz file on GENOME using HISAT2"
#	@echo "Example local usage"
#	@echo '    docker run -v "$(pwd):/data" -v "test/ref:/ref" unigebsp/rnaseq_pipeline GENOME=REF fqfile.star'
#	@echo "Example usage on Baobab:"
#	@echo "    # ask for a 16 cpu node for 8 hours"
#	@echo "    salloc -p shared-cpu --time=8:00:00 --mem=64G -n1 --cpus-per-task=16"
#	@echo "    # run the pipeline"
#	@echo "    singularity run docker://unigebsp/rnaseq-pipeline"


#-#-#-#-#-#-#-#-#
# Main rule
#-#-#-#-#-#-#-#-#
%/ALL:
	$(MAKE) "$*"/FASTQC "$*"/HT2

#-#-#-#-#-#-#-#-#
# testing
#-#-#-#-#-#-#-#-#
test:
	$(MAKE) GENOME=DdMm data/fastq/test/ALL
	$(MAKE) GENOME=Mm data/fastq/test/ALL
	$(MAKE) GENOME=Dd data/fastq/test/ALL
	

%/CLEAN:
	rm -f "$*"/*_fastqc.html "$*"/*.ht2.*

#-#-#-#-#-#-#-#-#-#-#-#-#-#
# import external rules
#-#-#-#-#-#-#-#-#-#-#-#-#-#
include $(LOCALDIR)/Makefile.samtools
include $(LOCALDIR)/Makefile.ht2
include $(LOCALDIR)/Makefile.fastq



#cwd := $(notdir $(patsubst %/,%,$(dir $(mkfile_path))))



