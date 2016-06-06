# dunetpc/dune/DataPrep/Service/test/CMakeLists.txt

# David Adams
# May 2016

# Build test for each service.

include(CetTest)

cet_test(test_StandardRawDigitExtractService SOURCES test_StandardRawDigitExtractService.cxx
  LIBRARIES
    dune_DataPrep_Service_StandardRawDigitExtractService_service
    dune_ArtSupport
)

cet_test(test_StandardRawDigitPrepService SOURCES test_StandardRawDigitPrepService.cxx
  LIBRARIES
    dune_DataPrep_Service_StandardRawDigitExtractService_service
    dune_ArtSupport
)



