# dunetpc/dune/DuneInterface/Data/test/CMakeLists.txt

# David Adams
# May 2018

# Build test for each component.

include(CetTest)

cet_test(test_DataMap SOURCES test_DataMap.cxx
  LIBRARIES
    ROOT_BASIC_LIB_LIST
)

cet_test(test_RunData SOURCES test_RunData.cxx
  LIBRARIES
    ROOT_BASIC_LIB_LIST
)

cet_test(test_AdcChannelData SOURCES test_AdcChannelData.cxx
  LIBRARIES
    ROOT_BASIC_LIB_LIST
)

cet_test(test_FftwReal2dDftData SOURCES test_FftwReal2dDftData.cxx
  LIBRARIES
    ROOT_BASIC_LIB_LIST
)

cet_test(test_Tpc2dRoi SOURCES test_Tpc2dRoi.cxx
  LIBRARIES
    ROOT_BASIC_LIB_LIST
)

cet_test(test_TpcData SOURCES test_TpcData.cxx
  LIBRARIES
    dune_DuneInterface_Data
    ROOT_BASIC_LIB_LIST
)

cet_enable_asserts()
