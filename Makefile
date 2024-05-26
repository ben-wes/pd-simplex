lib.name = simplex~
class.sources = simplex~.c
datafiles = simplex~-help.pd LICENSE
datadirs = x

PDLIBBUILDER_DIR=./pd-lib-builder
include $(PDLIBBUILDER_DIR)/Makefile.pdlibbuilder

# expecting pblibbuilder's $(system) definition
ifeq ($(system), Darwin)
  override c.ldflags = -undefined dynamic_lookup -flat_namespace -bundle
  override cxx.ldflags = -undefined dynamic_lookup -flat_namespace -bundle
endif
