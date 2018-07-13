# dunetpc/dune/Protodune/Tool/test/CMakeLists.txt

# David Adams
# July 2018

# Build test for each tool.

include(CetTest)

cet_test(test_ProtoDuneChannelRanges SOURCES test_ProtoDuneChannelRanges.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    ${ROOT_BASIC_LIB_LIST}
    ${SQLITE3}
    ${Boost_FILESYSTEM_LIBRARY}
)
