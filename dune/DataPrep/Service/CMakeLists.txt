# dunetpc/dune/DataPrep/Service/CMakeLists.txt

# David Adams
# May 2016

# Add plugin for each service.

simple_plugin(InterpolatingAdcMitigationService   "service"
                art_Utilities canvas_Utilities
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(StandardRawDigitExtractService   "service"
                lardataobj_RawData
                art_Utilities canvas_Utilities
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(ThresholdNoiseRemovalService   "service"
                art_Utilities canvas_Utilities
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(MultiChannelNoiseRemovalService   "service"
                art_Utilities canvas_Utilities
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(Dune35tNoiseRemovalService   "service"
                larcore_Geometry
                lbne-raw-data_Services_ChannelMap_ChannelMapService_service
                art_Utilities canvas_Utilities
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(MedianPedestalService   "service"
                art_Utilities canvas_Utilities
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcSuppressSignalFindingService "service"
                art_Utilities canvas_Utilities
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(DuneDeconvolutionService "service"
                art_Utilities canvas_Utilities
                ${ROOT_BASIC_LIB_LIST}
                SignalShapingServiceDUNE_service
             )

simple_plugin(DuneRoiBuildingService "service"
                art_Utilities canvas_Utilities
                ${ROOT_BASIC_LIB_LIST}
                SignalShapingServiceDUNE_service
             )

simple_plugin(StandardAdcWireBuildingService "service"
                lardataobj_RecoBase
                lardata_RecoBaseArt
                art_Utilities canvas_Utilities
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(StandardRawDigitPrepService   "service"
                lardataobj_RawData
                lardataobj_RecoBase
                art_Utilities canvas_Utilities
                ${ROOT_BASIC_LIB_LIST}
             )

add_subdirectory(test)
