# dunetpc/dune/DataPrep/Service/test/CMakeLists.txt

# David Adams
# May 2016

# Build test for each service.

include(CetTest)

cet_test(test_StandardRawDigitExtractService SOURCES test_StandardRawDigitExtractService.cxx
  LIBRARIES
    dune_DataPrep_Service_StandardRawDigitExtractService_service
    dune_ArtSupport
    ${ROOT_CORE}
)

cet_test(test_InterpolatingAdcMitigationService SOURCES test_InterpolatingAdcMitigationService.cxx
  LIBRARIES
    dune_ArtSupport
    ${ROOT_CORE}
)

cet_test(test_ThresholdNoiseRemovalService SOURCES test_ThresholdNoiseRemovalService
  LIBRARIES
    dune_ArtSupport
    ${ROOT_CORE}
)

cet_test(test_MultiChannelNoiseRemovalService SOURCES test_MultiChannelNoiseRemovalService
  LIBRARIES
    dune_ArtSupport
    ${ROOT_CORE}
)

cet_test(test_Dune35tNoiseRemovalService SOURCES test_Dune35tNoiseRemovalService
  LIBRARIES
    dune_ArtSupport
    lbne-raw-data_Services_ChannelMap_ChannelMapService_service
    ${ROOT_CORE}
)

cet_test(test_MedianPedestalService SOURCES test_MedianPedestalService
  LIBRARIES
    dune_DataPrep_Service_MedianPedestalService_service
    dune_ArtSupport
    ${ROOT_CORE}
)

cet_test(test_AdcSuppressSignalFindingService SOURCES test_AdcSuppressSignalFindingService
  LIBRARIES
    dune_DataPrep_Service_AdcSuppressSignalFindingService_service
    dune_ArtSupport
    ${ROOT_CORE}
)

cet_test(test_DuneDeconvolutionService SOURCES test_DuneDeconvolutionService
  LIBRARIES
    dune_DataPrep_Service_DuneDeconvolutionService_service
    dune_ArtSupport
    ${ROOT_CORE}
)

cet_test(test_DuneRoiBuildingService SOURCES test_DuneRoiBuildingService
  LIBRARIES
    dune_DataPrep_Service_DuneRoiBuildingService_service
    dune_ArtSupport
    ${ROOT_CORE}
)

cet_test(test_StandardAdcWireBuildingService SOURCES test_StandardAdcWireBuildingService
  LIBRARIES
    dune_DataPrep_Service_StandardAdcWireBuildingService_service
    dune_ArtSupport
    lardataobj_RawData
    ${ROOT_CORE}
)

cet_test(test_StandardRawDigitPrepService SOURCES test_StandardRawDigitPrepService.cxx
  LIBRARIES
    dune_DataPrep_Service_StandardRawDigitExtractService_service
    dune_DataPrep_Service_MedianPedestalService_service
    dune_ArtSupport
    ${ROOT_CORE}
)



