# Compile the utility programs in this package.

include SoftRelTools/arch_spec_root.mk
include SoftRelTools/arch_spec_pgsql.mk
include SoftRelTools/arch_spec_nutools.mk

LIB_TYPE        := shared
LIB             := lib$(PACKAGE)
LIBCXXFILES     := $(wildcard *.cxx)
JOBFILES        := $(wildcard *.fcl)

LIBLINK    := -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -l$(PACKAGE)  -L${PGSQL_DIR}/lib -lpq

########################################################################
include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_art.mk
########################################################################

override LIBLIBS += $(LOADLIBES) -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -L$(NUTOOLS_LIB) -lGeometry -lSimulationBase
