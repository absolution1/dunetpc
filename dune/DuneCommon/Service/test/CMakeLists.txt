# dunetpc/dune/DuneCommon/Service/test/CMakeLists.txt

# David Adams
# October 2016

# Build test for each service.

include(CetTest)

cet_test(test_GeoRopChannelGroupService SOURCES test_GeoRopChannelGroupService.cxx
  LIBRARIES
    dune_ArtSupport
)

cet_test(test_GeoApaChannelGroupService SOURCES test_GeoApaChannelGroupService.cxx
  LIBRARIES
    dune_ArtSupport
)
