
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Parameters
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
NJOB = 3
NCPU = 6
GENOMEDIR = ./data/ref
GENOME = Dd+Mm


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Set environment variables
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
MKDIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
export LOCALDIR ?= $(MKDIR)/.local/
export PATH := $(LOCALDIR)/bin/:$(PATH)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Main
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
.PHONY:usage %/all %/clean test

usage:
	@echo "usage: rnaseq <option> target"
	@echo "options:" 
	@echo "    -k                Continue on error"
	@echo "    -j <int>          Parallel execution on <int> jobs"
	@echo "    -n                Dry run: don't run anything, just show commands"
	@echo "    GENOMEDIR='path'  Directory containing the reference genomes [$(GENOMEDIR)]"
	@echo "    GENOME='path'     Default reference genome [$(GENOME)]"
	@echo "targets:"
	@echo "    <dir>/all               Run <dir>/FASTQC and <dir>/HT2"
	@echo "    <dir>/FASTQC            Run fastqc on all .fastq.gz files in the directory (recursively)"
	@echo "    <dir>/HT2               Map all .fastq.gz in the directory on the default genome (recursively)"	
	@echo "    fqfile.HT2              Map a single .fastq.gz file on GENOME using HISAT2"

%/all:
	$(MAKE) "$*"/FASTQC "$*"/HT2

%/clean:
	rm -rf "$*"/*_fastqc.html "$*"/*.ht2.*

test:
	$(MAKE) GENOME=Mm data/fastq/test/all
	$(MAKE) GENOME=Dd+Mm data/fastq/test/all
	$(MAKE) GENOME=Dd data/fastq/test/all
	


#-#-#-#-#-#-#-#-#-#-#-#-#-#
# import external rules
#-#-#-#-#-#-#-#-#-#-#-#-#-#
include $(LOCALDIR)/Makefile.ht2
include $(LOCALDIR)/Makefile.samtools
include $(LOCALDIR)/Makefile.fastq
include $(LOCALDIR)/Makefile.ref



#cwd := $(notdir $(patsubst %/,%,$(dir $(mkfile_path))))



