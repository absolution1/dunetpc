# dunetpc/dune/DuneCommon/Tool/test/CMakeLists.txt

# David Adams
# August 2017

# Build test for each tool.

include(CetTest)

cet_test(test_FclRunDataTool SOURCES test_FclRunDataTool.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    ROOT_BASIC_LIB_LIST
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
)

cet_test(test_FclIndexRangeTool SOURCES test_FclIndexRangeTool.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    ROOT_BASIC_LIB_LIST
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
)

cet_test(test_FclIndexVectorMapTool SOURCES test_FclIndexVectorMapTool.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    ROOT_BASIC_LIB_LIST
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
)

cet_test(test_FixedTimeOffsetTool SOURCES test_FixedTimeOffsetTool.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    ROOT_BASIC_LIB_LIST
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
)

cet_test(test_TimingRawDecoderOffsetTool SOURCES test_TimingRawDecoderOffsetTool.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    ROOT_BASIC_LIB_LIST
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
)

cet_test(test_SimpleHistogramManager SOURCES test_SimpleHistogramManager.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    ROOT_BASIC_LIB_LIST
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
)

cet_test(test_FclFloatArray SOURCES test_FclFloatArray.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    ROOT_BASIC_LIB_LIST
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
)
cet_test(test_FclFileFloatArray SOURCES test_FclFileFloatArray.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    ROOT_BASIC_LIB_LIST
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
)

cet_test(test_ChannelStatusServiceTool SOURCES test_ChannelStatusServiceTool.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    ROOT_BASIC_LIB_LIST
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
)

cet_test(test_ChannelStatusConfigTool SOURCES test_ChannelStatusConfigTool.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    ROOT_BASIC_LIB_LIST
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
)
cet_test(test_BlockIndexMapTool SOURCES test_BlockIndexMapTool.cxx
  LIBRARIES
    dune_ArtSupport
    art_Utilities
    canvas
    fhiclcpp
    cetlib cetlib_except
    cetlib_except
    ROOT_BASIC_LIB_LIST
    SQLITE3
    ${Boost_FILESYSTEM_LIBRARY}
)

