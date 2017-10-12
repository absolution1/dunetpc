# dunetpc/dune/DataPrep/Tool/test/CMakeLists.txt

# David Adams
# June 2017

# Build test for each service.

include(CetTest)

cet_test(test_AdcChannelDumper SOURCES test_AdcChannelDumper.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    ${SQLITE3}
    ${Boost_FILESYSTEM_LIBRARY}
)

cet_test(test_AdcChannelPlotter SOURCES test_AdcChannelPlotter.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    ${SQLITE3}
    ${Boost_FILESYSTEM_LIBRARY}
)

cet_test(test_AdcPedestalFitter SOURCES test_AdcPedestalFitter.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    ${SQLITE3}
    ${Boost_FILESYSTEM_LIBRARY}
)

cet_test(test_AdcDataDumper SOURCES test_AdcDataDumper.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    ${SQLITE3}
    ${Boost_FILESYSTEM_LIBRARY}
)

cet_test(test_AdcDataPlotter SOURCES test_AdcDataPlotter.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    ${SQLITE3}
    ${Boost_FILESYSTEM_LIBRARY}
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
    ${SQLITE3}
    ${Boost_FILESYSTEM_LIBRARY}
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
    ${SQLITE3}
    ${Boost_FILESYSTEM_LIBRARY}
)

cet_test(test_AdcThresholdSignalFinderAcdWireReader
  SOURCES
    test_AdcThresholdSignalFinder.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    ${SQLITE3}
    ${Boost_FILESYSTEM_LIBRARY}
)

cet_test(test_AcdWireReader SOURCES test_AcdWireReader.cxx
  LIBRARIES
    lardataobj_RawData
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    ${SQLITE3}
    ${Boost_FILESYSTEM_LIBRARY}
)

cet_enable_asserts()
