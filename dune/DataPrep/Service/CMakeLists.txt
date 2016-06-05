# dunetpc/dune/DataPrep/Service/CMakeLists.txt

# David Adams
# May 2016

# Add plugin for each service.

simple_plugin(StandardRawDigitExtractService   "service"
                lardataobj_RawData
                art_Utilities canvas_Utilities
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(StandardRawDigitPrepService   "service"
                lardataobj_RawData
                art_Utilities canvas_Utilities
                ${ROOT_BASIC_LIB_LIST}
             )

add_subdirectory(test)
