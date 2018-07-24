# $Id: GNUmakefile 42 2010-01-15 15:48:21Z vnivanch $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := tpcsim
G4TARGET := $(name)
G4EXLIB := true
###$(info Hello1)
.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk

#Add ROOT options for compilation
CPPFLAGS += `root-config --cflags`
LDFLAGS  += `root-config --libs ` -lGenVector

include $(G4INSTALL)/config/binmake.gmk

