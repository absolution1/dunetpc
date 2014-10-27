#
# $Id: GNUmakefile,v 1.2 2011/02/17 01:45:48 brebel Exp $
#
include SoftRelTools/arch_spec_root.mk

PACKAGE     := GeometryTest
LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
JOBFILES    := $(wildcard *.fcl)

LIBLINK     := -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(LARSOFT)/lib/$(SRT_SUBDIR) -l$(PACKAGE)


########################################################################
include SoftRelTools/standard.mk

include SoftRelTools/arch_spec_art.mk

override LIBLIBS += $(LOADLIBES) -L$(ROOTSYS)/lib -lHist -lRIO -lGeom -lXMLIO -lGdml -lTree -lGraf -lGui -lGed -lGeomBuilder -lGraf3d -lGpad -lGeomPainter -L$(ROOTSYS)/lib -lGeom -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(LARSOFT)/lib/$(SRT_SUBDIR) -lGeometry
