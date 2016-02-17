# dunetpc/dune/Daq/Service/test/CMakeLists.txt

include(CetTest)

MESSAGE( STATUS "CMAKE_LIBRARY_PATH: " ${CMAKE_LIBRARY_PATH} )

cet_test(test_FileChannelMappingService SOURCES test_FileChannelMappingService.cxx
  LIBRARIES
    dune_Daq_Service_FileChannelMappingService_service
    dune_ArtSupport
)
