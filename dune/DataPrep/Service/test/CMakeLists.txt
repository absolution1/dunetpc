# dunetpc/dune/DataPrep/Service/test/CMakeLists.txt

# David Adams
# May 2016

# Build test for each service.

include(CetTest)

cet_test(test_StandardRawDigitExtractService SOURCES test_StandardRawDigitExtractService.cxx
  LIBRARIES
    dune_DataPrep_Service_StandardRawDigitExtractService_service
    dune_ArtSupport
)

cet_test(test_InterpolatingAdcMitigationService SOURCES test_InterpolatingAdcMitigationService.cxx
  LIBRARIES
    dune_ArtSupport
)

cet_test(test_Dune35tNoiseRemovalService SOURCES test_Dune35tNoiseRemovalService
  LIBRARIES
    dune_ArtSupport
    lbne-raw-data_Services_ChannelMap_ChannelMapService_service
)

cet_test(test_MedianPedestalService SOURCES test_MedianPedestalService
  LIBRARIES
    dune_DataPrep_Service_MedianPedestalService_service
    dune_ArtSupport
)

cet_test(test_AdcSuppressSignalFindingService SOURCES test_AdcSuppressSignalFindingService
  LIBRARIES
    dune_DataPrep_Service_AdcSuppressSignalFindingService_service
    dune_ArtSupport
)

cet_test(test_DuneDeconvolutionService SOURCES test_DuneDeconvolutionService
  LIBRARIES
    dune_DataPrep_Service_DuneDeconvolutionService_service
    dune_ArtSupport
)

cet_test(test_DuneRoiBuildingService SOURCES test_DuneRoiBuildingService
  LIBRARIES
    dune_DataPrep_Service_DuneRoiBuildingService_service
    dune_ArtSupport
)

cet_test(test_StandardRawDigitPrepService SOURCES test_StandardRawDigitPrepService.cxx
  LIBRARIES
    dune_DataPrep_Service_StandardRawDigitExtractService_service
    dune_DataPrep_Service_MedianPedestalService_service
    dune_ArtSupport
)



