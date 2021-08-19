# dunetpc/dune/DataPrep/Tool/test/CMakeLists.txt

# David Adams
# June 2017

# Build test for each service.

include(CetTest)

cet_test(test_StandardAdcChannelStringTool SOURCES test_StandardAdcChannelStringTool.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_AdcResultDumper SOURCES test_AdcResultDumper.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_AdcChannelDumper SOURCES test_AdcChannelDumper.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_AdcChannelPlotter SOURCES test_AdcChannelPlotter.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_AdcTickModViewer SOURCES test_AdcTickModViewer.cxx
  LIBRARIES
    dune_DataPrep_Utility
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_AdcChannelTrimmer SOURCES test_AdcChannelTrimmer.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_AdcChannelSplitter SOURCES test_AdcChannelSplitter.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_AdcPedestalFitter SOURCES test_AdcPedestalFitter.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_AdcCodeMitigator SOURCES test_AdcCodeMitigator.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_ExpTailRemover SOURCES test_ExpTailRemover.cxx
  LIBRARIES
    dune_ArtSupport
    dune_DuneCommon_Utility
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_ExpTailPedRemover SOURCES test_ExpTailPedRemover.cxx
  LIBRARIES
    dune_ArtSupport
    dune_DuneCommon_Utility
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_VintageDeconvoluter SOURCES test_VintageDeconvoluter.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_AdcSampleScaler SOURCES test_AdcSampleScaler.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_AdcRangeSampleScaler SOURCES test_AdcRangeSampleScaler.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_AdcDataDumper SOURCES test_AdcDataDumper.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_AdcEventViewer SOURCES test_AdcEventViewer.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_AdcDataPlotter SOURCES test_AdcDataPlotter.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_AdcDetectorPlotter SOURCES test_AdcDetectorPlotter.cxx
  LIBRARIES
    ART_FRAMEWORK_SERVICES_REGISTRY
    dune_ArtSupport
    dune_Geometry
    larcorealg_Geometry
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST

)

cet_test(test_AdcChannelMetric SOURCES test_AdcChannelMetric.cxx
  LIBRARIES
    dune_ArtSupport
    dune_Geometry
    larcorealg_Geometry
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_AcdDigitReader SOURCES test_AcdDigitReader.cxx
  LIBRARIES
    lardataobj_RawData
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_AdcSampleFiller
  SOURCES
    test_AdcSampleFiller.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_FclStickyCodeFlagger
  SOURCES
    test_FclStickyCodeFlagger.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_AdcKeepAllSignalFinder
  SOURCES
    test_AdcKeepAllSignalFinder.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_AdcThresholdSignalFinder
  SOURCES
    test_AdcThresholdSignalFinder.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_AdcNoiseSignalFinder
  SOURCES
    test_AdcNoiseSignalFinder.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_AdcRoiShifter
  SOURCES
    test_AdcRoiShifter.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_AdcRegularSignalFinder
  SOURCES
    test_AdcRegularSignalFinder.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_DuneAdcSignalFinder
  SOURCES
    test_DuneAdcSignalFinder.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_AcdWireReader SOURCES test_AcdWireReader.cxx
  LIBRARIES
    ART_FRAMEWORK_SERVICES_REGISTRY
    lardataobj_RawData
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_AdcRoiViewer
  SOURCES
    test_AdcRoiViewer.cxx
  LIBRARIES
    lardataobj_RawData
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    ROOT_BASIC_LIB_LIST
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
)

cet_test(test_FembLinearCalibration
  SOURCES
    test_FembLinearCalibration.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_FloatArrayGainCalibration
  SOURCES
    test_FloatArrayGainCalibration
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_AdcChannelFFT
  SOURCES
    test_AdcChannelFFT.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_Adc2dConvolute SOURCES test_Adc2dConvolute.cxx
  LIBRARIES
    dune_ArtSupport
    dune_DuneCommon_Utility
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_AdcDeconvoluteFFT
  SOURCES
    test_AdcDeconvoluteFFT.cxx
  LIBRARIES
    dune_ArtSupport
    dune_DuneCommon_Utility
    art_Utilities
    canvas
    fhiclcpp
    ROOT_BASIC_LIB_LIST
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
)

cet_test(test_AdcChannelDftPlotter
  SOURCES
    test_AdcChannelDftPlotter.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_FembMappedAdcModifier
  SOURCES
    test_FembMappedAdcModifier.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_AdcRoiSlicer
  SOURCES
    test_AdcRoiSlicer.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_test(test_AdcRoiToTree
  SOURCES
    test_AdcRoiToTree.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

cet_enable_asserts()
