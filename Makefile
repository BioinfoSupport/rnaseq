
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Parameters
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
NJOB = 3
NCPU = 6
GENOMEDIR = ./data/ref
GENOME = Mm


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Set environment variables
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
MKDIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
export LOCALDIR ?= $(MKDIR)/.local/
export PATH := $(LOCALDIR)/bin/:$(PATH)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Main
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
.PHONY:USAGE %/ALL %/CLEAN TEST

USAGE:
	@echo "usage: rnaseq <option> target"
	@echo "options:" 
	@echo "    -k                Continue on error"
	@echo "    -j <int>          Parallel execution on <int> jobs"
	@echo "    -n                Dry run: don't run anything, just show commands"
	@echo "    GENOMEDIR='path'  Directory containing the reference genomes [$(GENOMEDIR)]"
	@echo "    GENOME='path'     Default reference genome [$(GENOME)]"
	@echo "targets:"




TEST:
	$(MAKE) GENOME=Mm data/fastq/test/RNASEQ
	#$(MAKE) GENOME=Dd+Mm data/fastq/test/ALL
	#$(MAKE) GENOME=Dd data/fastq/test/ALL
	

%/CLEAN:
	rm -rf $*/*.ht2.* $*/*.fastq.gz.* $*/*.html

#-#-#-#-#-#-#-#-#-#-#-#-#-#
# import external rules
#-#-#-#-#-#-#-#-#-#-#-#-#-#
include $(LOCALDIR)/Makefile.rnaseq



#cwd := $(notdir $(patsubst %/,%,$(dir $(mkfile_path))))



