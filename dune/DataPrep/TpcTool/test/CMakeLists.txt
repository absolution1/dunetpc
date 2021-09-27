# dunetpc/dune/DataPrep/TpcTool/test/CMakeLists.txt

# David Adams
# February 2021

# Build test for each tool.

include(CetTest)

cet_test(test_AdcToRoi2d SOURCES test_AdcToRoi2d.cxx
  LIBRARIES
    dune_DuneInterface_Data
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

cet_test(test_Roi2dToAdc SOURCES test_Roi2dToAdc.cxx
  LIBRARIES
    dune_DuneInterface_Data
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

cet_test(test_Tpc2dDeconvolute SOURCES test_Tpc2dDeconvolute.cxx
  LIBRARIES
    dune_DuneCommon_Utility
    dune_DuneInterface_Data
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    ${SQLITE3}
    ${Boost_FILESYSTEM_LIBRARY}
    ${ROOT_BASIC_LIB_LIST}
)

cet_enable_asserts()
