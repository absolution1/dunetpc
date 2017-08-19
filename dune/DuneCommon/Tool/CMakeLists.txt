# dunetpc/dune/DuneCommon/Tool/CMakeLists.txt

# David Adams
# August 2017

# Add plugin for each tool.

simple_plugin(SimpleHistogramManager  "tool"
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

add_subdirectory(test)
