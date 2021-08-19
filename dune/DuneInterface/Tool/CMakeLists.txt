# dunetpc/dune/DuneInterface/ToolCMakeLists.txt
#
# David Adams
# August 2017
#
# Chris Green
# August 2021
#
# CMake instructions for DUNE tool interface classes.

# Make the dictionary.
art_dictionary()


####################################
# Define interface targets (no code, only depedencies) when headers do
# not have an associated implementation compilation unit, but do confer
# a transitive dependency on other targets.
cet_make_library(LIBRARY_NAME AdcChannelTool INTERFACE
  SOURCE AdcChannelTool.h
  LIBRARIES INTERFACE dune::DuneInterface_Data_DataMap)

cet_make_library(LIBRARY_NAME AdcChannelStringTool INTERFACE
  SOURCE AdcChannelStringTool.h
  LIBRARIES INTERFACE
  dune::DuneInterface_Data
  )

cet_make_library(LIBRARY_NAME RunDataTool INTERFACE
  SOURCE RunDataTool.h
  LIBRARIES INTERFACE
  dune::DuneInterface_Data_RunData
  )

cet_make_library(LIBRARY_NAME TpcDataTool INTERFACE
  SOURCE TpcDataTool.h
  LIBRARIES INTERFACE
  dune::DuneInterface_Data
  dune::AdcChannelTool
)


install_headers()
install_source( README.txt )

add_subdirectory(test)
