# dunetpc/dune/ArtSupport/Tool/CMakeLists.txt

# David Adams
# August 2017

# Add plugin for each tool.

simple_plugin(TestTool   "tool"
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )
