# dunetpc/dune/DataPrep/Tool/CMakeLists.txt

# David Adams
# June 2017

# Add plugin for each tool.

simple_plugin(AdcDataDumper   "tool"
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcDataPlotter   "tool"
                dune_DuneCommon
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AcdDigitReader   "tool"
                dune_DuneCommon
                lardataobj_RawData
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

add_subdirectory(test)
