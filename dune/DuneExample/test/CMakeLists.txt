# dunetpc/dune/DuneCommon/test/CMakeLists.txt

include(CetTest)

cet_enable_asserts()

cet_test(test_MyClass SOURCES test_MyClass.cxx
  LIBRARIES
    dune_DuneExample
    ROOT_BASIC_LIB_LIST
)

