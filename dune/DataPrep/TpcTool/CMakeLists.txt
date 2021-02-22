# dunetpc/dune/DataPrep/TpcTool/CMakeLists.txt

# David Adams
# February 2021
#
# Tools that act on TpcData.
#
# Add plugin for each tool.

simple_plugin(AdcToRoi2d   "tool"
                dune_DuneInterface_Data
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(Roi2dToAdc   "tool"
                dune_DuneInterface_Data
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

add_subdirectory(test)
install_source()

