# dunetpc/dune/DuneInterface/Data/test/CMakeLists.txt

# David Adams
# May 2018

# Build test for each component.

include(CetTest)

cet_test(test_DataMap SOURCES test_DataMap.cxx
  LIBRARIES
    ${ROOT_BASIC_LIB_LIST}
)

cet_test(test_RunData SOURCES test_RunData.cxx
  LIBRARIES
    ${ROOT_BASIC_LIB_LIST}
)

cet_enable_asserts()
