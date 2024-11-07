lib.name = simplex
class.sources = simplex~.c simplex.c
datafiles = simplex~-help.pd simplex-help.pd simplex~-meta.pd LICENSE
datadirs = x

PDLIBBUILDER_DIR=./pd-lib-builder
include $(PDLIBBUILDER_DIR)/Makefile.pdlibbuilder
