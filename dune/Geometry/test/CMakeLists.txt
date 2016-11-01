# dunetpc/dune/Geometry/test/CMakeLists.txt

# David Adams
# October 2016

# Build tests for channel maps.

include(CetTest)

cet_test(test_Geometry_Dune35t SOURCES test_Geometry_Dune35t.cxx
  LIBRARIES
    dune_ArtSupport
    larcore_Geometry
)
