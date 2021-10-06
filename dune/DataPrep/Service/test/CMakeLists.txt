# dunetpc/dune/DataPrep/Service/test/CMakeLists.txt

# David Adams
# May 2016

# Build test for each service.

include(CetTest)

cet_enable_asserts()

cet_test(test_ConfigurableAdcChannelDataCopyService SOURCES test_ConfigurableAdcChannelDataCopyService.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    ART_FRAMEWORK_SERVICES_REGISTRY
    cetlib
    cetlib_except
    ROOT_BASIC_LIB_LIST
)

cet_test(test_StandardRawDigitExtractService SOURCES test_StandardRawDigitExtractService.cxx
  LIBRARIES
    dune_ArtSupport
    lardataobj_RawData
    art_Utilities
    canvas
    ART_FRAMEWORK_SERVICES_REGISTRY
    cetlib
    cetlib_except
    ROOT_BASIC_LIB_LIST
)

cet_test(test_InterpolatingAdcMitigationService SOURCES test_InterpolatingAdcMitigationService.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    ART_FRAMEWORK_SERVICES_REGISTRY
    cetlib
    cetlib_except
    ROOT_BASIC_LIB_LIST
)

cet_test(test_ThresholdNoiseRemovalService SOURCES test_ThresholdNoiseRemovalService
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    ART_FRAMEWORK_SERVICES_REGISTRY
    cetlib
    cetlib_except
    ROOT_BASIC_LIB_LIST
)

cet_test(test_MultiChannelNoiseRemovalService SOURCES test_MultiChannelNoiseRemovalService
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    ART_FRAMEWORK_SERVICES_REGISTRY
    cetlib
    cetlib_except
    ROOT_BASIC_LIB_LIST
)

cet_test(test_MedianPedestalService SOURCES test_MedianPedestalService
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    cetlib
    cetlib_except
    ART_FRAMEWORK_SERVICES_REGISTRY
    cetlib
    cetlib_except
    ROOT_BASIC_LIB_LIST
)

cet_test(test_AdcSuppressSignalFindingService SOURCES test_AdcSuppressSignalFindingService
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    ART_FRAMEWORK_SERVICES_REGISTRY
    cetlib
    cetlib_except
    ROOT_BASIC_LIB_LIST
)

cet_test(test_DuneDeconvolutionService SOURCES test_DuneDeconvolutionService
  LIBRARIES
    SignalShapingServiceDUNE_service
    dune_ArtSupport
    art_Utilities
    canvas
    ART_FRAMEWORK_SERVICES_REGISTRY
    cetlib
    cetlib_except
    ROOT_BASIC_LIB_LIST
)

cet_test(test_KeepAllRoiBuildingService SOURCES test_KeepAllRoiBuildingService
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    ART_FRAMEWORK_SERVICES_REGISTRY
    cetlib
    cetlib_except
    ROOT_BASIC_LIB_LIST
)

cet_test(test_DuneRoiBuildingService SOURCES test_DuneRoiBuildingService
  LIBRARIES
    dune_DataPrep_Service_DuneRoiBuildingService_service
    dune_ArtSupport
    art_Utilities
    canvas
    ART_FRAMEWORK_SERVICES_REGISTRY
    cetlib
    cetlib_except
    ROOT_BASIC_LIB_LIST
)

cet_test(test_StandardAdcWireBuildingService SOURCES test_StandardAdcWireBuildingService
  LIBRARIES
    dune_ArtSupport
    lardataobj_RawData
    art_Utilities
    canvas
    ART_FRAMEWORK_SERVICES_REGISTRY
    cetlib
    cetlib_except
    ROOT_BASIC_LIB_LIST
)

cet_test(test_StandardRawDigitPrepService SOURCES test_StandardRawDigitPrepService.cxx
  LIBRARIES
    dune_ArtSupport
    lardataobj_RawData
    art_Utilities
    canvas
    ART_FRAMEWORK_SERVICES_REGISTRY
    cetlib
    cetlib_except
    ROOT_BASIC_LIB_LIST
)

cet_test(test_ToolBasedRawDigitPrepService SOURCES test_ToolBasedRawDigitPrepService.cxx
  LIBRARIES
    dune_ArtSupport
    lardataobj_RawData
    art_Utilities
    canvas
    ROOT_HIST
    ART_FRAMEWORK_SERVICES_REGISTRY
    cetlib
    cetlib_except
    ROOT_BASIC_LIB_LIST
)

cet_test(test_TpcToolBasedRawDigitPrepService SOURCES test_TpcToolBasedRawDigitPrepService.cxx
  LIBRARIES
    dune_ArtSupport
    lardataobj_RawData
    art_Utilities
    canvas
    ROOT_HIST
    ART_FRAMEWORK_SERVICES_REGISTRY
    cetlib
    cetlib_except
    ROOT_BASIC_LIB_LIST
)
