
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Set environment variables
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
MKDIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
export LOCALDIR ?= $(MKDIR)/.local/
export PATH := $(LOCALDIR)/bin/:$(PATH)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Parameters
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#include $(LOCALDIR)/config/Makefile.baobab
NJOB = 3
NCPU = 6
GENOMEDIR = ./data/ref
GENOME = Mm


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Main
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
.PHONY:USAGE
USAGE:
	@echo "usage: rnaseq <option> target"
	@echo "options:" 
	@echo "    -k                Continue on error"
	@echo "    -j <int>          Parallel execution on <int> jobs"
	@echo "    -n                Dry run: don't run anything, just show commands"
	@echo "    GENOMEDIR='path'  Directory containing the reference genomes [$(GENOMEDIR)]"
	@echo "    GENOME='path'     Default reference genome [$(GENOME)]"
	@echo "targets:"


.PHONY:TESTS
TESTS:
	$(MAKE) GENOME=Mm data/fastq/test/RNASEQ.ALL
	$(MAKE) GENOME=Dd+Mm data/fastq/test/RNASEQ.ALL
	$(MAKE) GENOME=Dd data/fastq/test/RNASEQ.ALL
	

.PHONY:%/CLEAN
%/CLEAN:
	rm -rf $*/*.ht2.* $*/*.fastq.gz.* $*/*.html

#-#-#-#-#-#-#-#-#-#-#-#-#-#
# import external rules
#-#-#-#-#-#-#-#-#-#-#-#-#-#
include $(LOCALDIR)/Makefile.rnaseq
#include $(LOCALDIR)/Makefile.tnseq



#cwd := $(notdir $(patsubst %/,%,$(dir $(mkfile_path))))



