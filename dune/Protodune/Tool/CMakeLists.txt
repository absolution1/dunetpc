# dunetpc/dune/Protodune/Tool/CMakeLists.txt

# David Adams
# July 2018

# Add plugin for each tool.

simple_plugin(ProtoDuneChannelRanges  "tool"
                dune_ArtSupport
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

add_subdirectory(test)
