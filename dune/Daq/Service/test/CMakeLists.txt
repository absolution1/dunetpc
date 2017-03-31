# dunetpc/dune/Daq/Service/test/CMakeLists.txt

include(CetTest)

cet_enable_asserts()

MESSAGE( STATUS "CMAKE_LIBRARY_PATH: " ${CMAKE_LIBRARY_PATH} )

# Comment test out b/c the channel map is not osolete.
# See https://cdcvs.fnal.gov/redmine/issues/12788
#cet_test(test_FileChannelMappingService SOURCES test_FileChannelMappingService.cxx
#  LIBRARIES
#    dune_Daq_Service_FileChannelMappingService_service
#    dune_ArtSupport
#)
