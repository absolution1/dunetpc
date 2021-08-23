# dunetpc/dune/Protodune/singlephase/Tool/CMakeLists.txt

# David Adams
# June 2018

# Add plugin for each tool.

cet_build_plugin(PdspOnlineChannel   art::tool
                dune_raw_data::dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
                art_Utilities canvas
                cetlib cetlib_except
                ART_FRAMEWORK_SERVICES_REGISTRY
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(IcebergOnlineChannel   art::tool
                dune_raw_data::dune-raw-data_Services_ChannelMap_IcebergChannelMapService_service
                art_Utilities canvas
                cetlib cetlib_except
                ART_FRAMEWORK_SERVICES_REGISTRY
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(ProtoduneOnlineChannel   art::tool
                art_Utilities canvas
                cetlib cetlib_except
                ART_FRAMEWORK_SERVICES_REGISTRY
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(ProtoDuneChannelRanges  art::tool
                dune_ArtSupport
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(ApaChannelRanges  art::tool
                dune_ArtSupport
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(HalfApaChannelRanges  art::tool
                dune_ArtSupport
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(ProtoDuneChannelGroups  art::tool
                dune_ArtSupport
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(ApaChannelGroups  art::tool
                dune_ArtSupport
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

add_subdirectory(test)
