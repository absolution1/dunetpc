# dunetpc/dune/DuneCommon/test/CMakeLists.txt

include(CetTest)

cet_enable_asserts()

cet_test(test_DuneTimeConverter SOURCES test_DuneTimeConverter.cxx
  LIBRARIES
    lbne_daqinput35t
    art_Persistency_Provenance canvas
)

cet_test(test_TPadManipulator SOURCES test_TPadManipulator.cxx
  LIBRARIES
    dune_DuneCommon
    ${ROOT_BASIC_LIB_LIST}
)
