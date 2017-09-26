# dunetpc/dune/DuneInterface/Data/test/CMakeLists.txt

# David Adams
# September 2017

# Build test for each component.

include(CetTest)

cet_test(test_DataMap SOURCES test_DataMap.cxx)

cet_enable_asserts()
