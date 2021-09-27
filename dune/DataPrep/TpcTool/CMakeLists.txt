# dunetpc/dune/DataPrep/TpcTool/CMakeLists.txt

# David Adams
# February 2021
#
# Tools that act on TpcData.
#
# Add plugin for each tool.

cet_build_plugin(AdcToRoi2d   art::tool
                dune_DuneInterface_Data
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(Roi2dToAdc   art::tool
                dune_DuneInterface_Data
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

simple_plugin(Tpc2dDeconvolute   "tool"
                dune_DuneInterface_Data
                dune_DuneCommon_Utility
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

add_subdirectory(test)
install_source()

