# dunetpc/dune/DuneCommon/test/CMakeLists.txt

include(CetTest)

cet_test(test_DuneTimeConverter SOURCES test_DuneTimeConverter.cxx
  LIBRARIES
    lbne_daqinput35t
    art_Persistency_Provenance canvas_Persistency_Provenance
)
