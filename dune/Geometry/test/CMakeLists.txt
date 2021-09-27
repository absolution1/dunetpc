# dunetpc/dune/Geometry/test/CMakeLists.txt

# David Adams
# November 2016

# Build tests for channel maps.
install_fhicl()
include(CetTest)
add_subdirectory(testIssue19191)

cet_enable_asserts()

cet_test(test_GeometryDune35t SOURCES test_GeometryDune35t.cxx
  LIBRARIES
    dune_ArtSupport
    larcorealg_Geometry
    ROOT_GEOM
    ROOT_BASIC_LIB_LIST
)

cet_test(test_GeometryDune10kt SOURCES test_GeometryDune10kt.cxx
  LIBRARIES
    dune_ArtSupport
    larcorealg_Geometry
    ROOT_GEOM
    ROOT_BASIC_LIB_LIST
)

cet_test(test_GeometryProtoDune SOURCES test_GeometryProtoDune.cxx
  LIBRARIES
    dune_ArtSupport
    larcorealg_Geometry
    ROOT_GEOM
    ROOT_BASIC_LIB_LIST
)

cet_test(test_WireSelector SOURCES test_WireSelector.cxx
  LIBRARIES
    dune_ArtSupport
    dune_Geometry
    dune_DuneCommon_Utility
    larcorealg_Geometry
    ROOT_GEOM
    ROOT_BASIC_LIB_LIST
)
