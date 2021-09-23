
include_directories( "$ENV{DUNEPDSPRCE_INC}" ) 
cet_find_library( RCEDAMLIB NAMES protodune-dam PATHS ENV DUNEPDSPRCE_LIB NO_DEFAULT_PATH )

cet_build_plugin(IcebergTPCRawDecoder art::module
                        lardataobj_RawData
                        dune_raw_data::dune-raw-data_Overlays
			artdaq_core::artdaq-core_Data
                        artdaq-core_Utilities
                        dune_raw_data::dune-raw-data_Services_ChannelMap_IcebergChannelMapService_service
                        ART_FRAMEWORK_CORE
                        ART_FRAMEWORK_PRINCIPAL
                        ART_FRAMEWORK_SERVICES_REGISTRY
                        ART_ROOT_IO_TFILE_SUPPORT
                        ROOT::Core
                        ART_ROOT_IO_TFILESERVICE_SERVICE
 #                       art_Persistency_Common canvas_Persistency_Common
 #                       art_Persistency_Provenance canvas_Persistency_Provenance
 #                       art_Utilities canvas_Utilities
                        ART_PERSISTENCY_PROVENANCE
                        messagefacility::MF_MessageLogger
                        CETLIB
                        ROOT_BASIC_LIB_LIST
                        RCEDAMLIB
			z
                        BASENAME_ONLY
)


cet_build_plugin(IcebergFELIXBufferDecoderMarch2021 art::module
                        lardataobj_RawData
                        dune_raw_data::dune-raw-data_Overlays
			artdaq_core::artdaq-core_Data
                        artdaq-core_Utilities
                        dune_raw_data::dune-raw-data_Services_ChannelMap_IcebergChannelMapService_service
                        ART_FRAMEWORK_CORE
                        ART_FRAMEWORK_PRINCIPAL
                        ART_FRAMEWORK_SERVICES_REGISTRY
                        ART_ROOT_IO_TFILE_SUPPORT
                        ROOT::Core
                        ART_ROOT_IO_TFILESERVICE_SERVICE
 #                       art_Persistency_Common canvas_Persistency_Common
 #                       art_Persistency_Provenance canvas_Persistency_Provenance
 #                       art_Utilities canvas_Utilities
                        ART_PERSISTENCY_PROVENANCE
                        messagefacility::MF_MessageLogger
                        CETLIB
                        ROOT_BASIC_LIB_LIST
                        RCEDAMLIB
			z
                        BASENAME_ONLY
)

cet_build_plugin(IcebergDataInterface   art::tool
				     art_utilities
				     canvas
				     cetlib
				     cetlib_except
                        lardataobj_RawData
                        dune_raw_data::dune-raw-data_Overlays
			artdaq_core::artdaq-core_Data
			artdaq-core_Utilities
                        dune_raw_data::dune-raw-data_Services_ChannelMap_IcebergChannelMapService_service
                        ART_FRAMEWORK_CORE
                        ART_FRAMEWORK_PRINCIPAL
                        ART_FRAMEWORK_SERVICES_REGISTRY
                        ART_ROOT_IO_TFILE_SUPPORT
                        ROOT::Core
                        ART_ROOT_IO_TFILESERVICE_SERVICE
                        ART_PERSISTENCY_PROVENANCE
                        messagefacility::MF_MessageLogger
                        CETLIB
                        ROOT_BASIC_LIB_LIST
                        RCEDAMLIB
                        z
             )


cet_build_plugin(IcebergDataInterfaceFELIXBufferMarch2021   art::tool
				     art_utilities
				     canvas
				     cetlib
				     cetlib_except
                        lardataobj_RawData
                        dune_raw_data::dune-raw-data_Overlays
			artdaq_core::artdaq-core_Data
			artdaq-core_Utilities
                        dune_raw_data::dune-raw-data_Services_ChannelMap_IcebergChannelMapService_service
                        ART_FRAMEWORK_CORE
                        ART_FRAMEWORK_PRINCIPAL
                        ART_FRAMEWORK_SERVICES_REGISTRY
                        ART_ROOT_IO_TFILE_SUPPORT
                        ROOT::Core
                        ART_ROOT_IO_TFILESERVICE_SERVICE
                        ART_PERSISTENCY_PROVENANCE
                        messagefacility::MF_MessageLogger
                        CETLIB
                        ROOT_BASIC_LIB_LIST
                        RCEDAMLIB
                        z
             )


install_headers()
install_fhicl()
install_source()
install_scripts()
