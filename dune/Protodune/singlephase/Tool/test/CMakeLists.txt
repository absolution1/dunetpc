# dunetpc/dune/Protodune/singlephase/Tool/test/CMakeLists.txt

# David Adams
# June 2018

# Build test for each tool.

include(CetTest)

cet_test(test_PdspOnlineChannel SOURCES test_PdspOnlineChannel.cxx
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

cet_test(test_ProtoduneOnlineChannel SOURCES test_ProtoduneOnlineChannel.cxx
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
