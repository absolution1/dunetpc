# art_extensions/test/AXService/CMakeLists.txt

include(CetTest)

cet_test(test_ArtServiceHelper SOURCES test_ArtServiceHelper.cxx
  LIBRARIES
    dune_ArtSupport
    art_Framework_Services_Optional_TFileService_service
    art_Framework_Services_Optional_RandomNumberGenerator_service
    ${ROOT_BASIC_LIB_LIST}
)

install_source()
