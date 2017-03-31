# dunetpc/dune/daqinput35t/test/CMakeLists.txt

include(CetTest)

cet_enable_asserts()

cet_test(test_tpcFragmentToRawDigits SOURCES test_tpcFragmentToRawDigits.cxx
  LIBRARIES
    lbne_daqinput35t
    art_Persistency_Provenance canvas
)
