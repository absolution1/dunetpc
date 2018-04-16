# dunetpc/dune/DuneCommon/Service/test/CMakeLists.txt

# David Adams
# October 2016

# Build test for each service.

include(CetTest)

cet_enable_asserts()

cet_test(test_GeoRopChannelGroupService
  SOURCES
    test_GeoRopChannelGroupService.cxx
  LIBRARIES
    larcorealg_Geometry
    dune_ArtSupport
)

cet_test(test_GeoApaChannelGroupService
  SOURCES
    test_GeoApaChannelGroupService.cxx
  LIBRARIES
    larcorealg_Geometry
    dune_ArtSupport
)

cet_test(test_FixedChannelGroupService
  SOURCES
    test_FixedChannelGroupService.cxx
  LIBRARIES
    dune_ArtSupport
    dune_DuneServiceAccess
)

cet_test(test_Dune35tChannelGroupService
  SOURCES
    test_Dune35tChannelGroupService.cxx
  LIBRARIES
    dune_ArtSupport
)

