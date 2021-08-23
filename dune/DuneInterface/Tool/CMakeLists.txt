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

###################################
# Generate CMake tool builder functions for use by build_plugin().

include(MakeDuneToolBuilder)
make_dune_tool_builder(dune::AdcChannelStringTool)
make_dune_tool_builder(dune::FloatArrayTool)
make_dune_tool_builder(dune::IndexMapTool)
make_dune_tool_builder(dune::IndexRangeGroupTool)
make_dune_tool_builder(dune::IndexRangeTool)
make_dune_tool_builder(dune::IndexVectorMapTool)
make_dune_tool_builder(dune::RunDataTool)
make_dune_tool_builder(dune::TimeOffsetTool)
make_dune_tool_builder(dune::TpcDataTool)

cet_collect_plugin_builders(Modules dune::ToolTypes)
####################################

install_headers()
install_source()

add_subdirectory(test)
