# dunetpc/dune/daqinput35t/test/CMakeLists.txt

include(CetTest)

cet_test(test_DuneTimeConverter SOURCES test_DuneTimeConverter.cxx
  LIBRARIES
    lbne_daqinput35t
    ${ART_PERSISTENCY_PROVENANCE}
)
