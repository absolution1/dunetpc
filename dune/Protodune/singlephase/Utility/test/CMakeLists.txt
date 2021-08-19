# dunetpc/dune/DataPrep/Utility/test/CMakeLists.txt

# David Adams
# June 2018

# Build test for each utility.

include(CetTest)

cet_test(test_ProtoduneChannelHelper SOURCES test_ProtoduneChannelHelper.cxx
  LIBRARIES
    dune_Protodune_singlephase_Utility
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
)

