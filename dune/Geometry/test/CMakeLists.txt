# dunetpc/dune/Geometry/test/CMakeLists.txt

# David Adams
# November 2016

# Build tests for channel maps.

include(CetTest)

cet_enable_asserts()

cet_test(test_GeometryDune35t SOURCES test_GeometryDune35t.cxx
  LIBRARIES
    dune_ArtSupport
    larcore_Geometry
)

cet_test(test_GeometryDune10kt SOURCES test_GeometryDune10kt.cxx
  LIBRARIES
    dune_ArtSupport
    larcore_Geometry
)

cet_test(test_GeometryProtoDune SOURCES test_GeometryProtoDune.cxx
  LIBRARIES
    dune_ArtSupport
    larcore_Geometry
)
