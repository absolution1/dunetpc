# dunetpc/dune/DataPrep/Service/CMakeLists.txt

# David Adams
# May 2016

# Add plugin for each service.

simple_plugin(InterpolatingAdcMitigationService   "service"
                ${ART_UTILITIES}
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(StandardRawDigitExtractService   "service"
                lardata_RawData
                ${ART_UTILITIES}
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(Dune35tNoiseRemovalService   "service"
                larcore_Geometry
                lbne-raw-data_Services_ChannelMap_ChannelMapService_service
                ${ART_UTILITIES}
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(MedianPedestalService   "service"
                ${ART_UTILITIES}
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(StandardRawDigitPrepService   "service"
                lardata_RawData
                ${ART_UTILITIES}
                ${ROOT_BASIC_LIB_LIST}
             )

add_subdirectory(test)
