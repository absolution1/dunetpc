# dunetpc/dune/DataPrep/Service/CMakeLists.txt

# David Adams
# May 2016

# Add plugin for each service.

simple_plugin(StandardRawDigitExtractService   "service"
                lardata_RawData
                ${ART_UTILITIES}
                ${ROOT_BASIC_LIB_LIST}
             )

add_subdirectory(test)
