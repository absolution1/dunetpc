# dunetpc/dune/DataPrep/Tool/test/CMakeLists.txt

# David Adams
# June 2017

# Build test for each service.

include(CetTest)

cet_enable_asserts()

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
