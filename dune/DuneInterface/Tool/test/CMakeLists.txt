# dunetpc/dune/DuneInterface/Tool/test/CMakeLists.txt

# David Adams
# April 2018

# Build test for each tool interface.

include(CetTest)

cet_test(test_AdcChannelTool SOURCES test_AdcChannelTool.cxx
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
