include_directories( "$ENV{DUNEPDSPRCE_INC}" )
cet_find_library( RCEDAMLIB NAMES protodune-dam PATHS ENV DUNEPDSPRCE_LIB NO_DEFAULT_PATH )

cet_build_plugin(PDSPTPCDataInterface   art::tool
				     art_utilities
				     canvas
				     cetlib
				     cetlib_except
                        lardataobj_RawData
                        dune_raw_data::dune-raw-data_Overlays
			artdaq_core::artdaq-core_Data
			artdaq-core_Utilities
                        dune_raw_data::dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
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


cet_build_plugin(TimingRawDecoder art::module
                        lardataobj_RawData
                        dune_raw_data::dune-raw-data_Overlays
                        dune_DuneObj
			artdaq-core_Utilities
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
                        BASENAME_ONLY
)

cet_build_plugin(SSPRawDecoder art::module
                        lardataobj_RawData
                        lardataobj_RecoBase
                        dune_raw_data::dune-raw-data_Overlays
                        dune_raw_data::dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
			artdaq_core::artdaq-core_Data
			artdaq-core_Utilities
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
                        BASENAME_ONLY
)

#[[
cet_build_plugin(PDWaveform art::module
                         lardataobj_RawData
                         lardataobj_RecoBase
                         dune_raw_data::dune-raw-data_Overlays
                         dune_raw_data::dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
 			artdaq_core::artdaq-core_Data
			artdaq-core_Utilities
                         ART_FRAMEWORK_CORE
                         ART_FRAMEWORK_PRINCIPAL
                         ART_FRAMEWORK_SERVICES_REGISTRY
                         ART_ROOT_IO_TFILE_SUPPORT
                         ROOT::Core
                         ART_ROOT_IO_TFILESERVICE_SERVICE
                         messagefacility::MF_MessageLogger
			 dune_Geometry
			 CETLIB
                         ROOT_BASIC_LIB_LIST
			 ROOT_SPECTRUM
                         BASENAME_ONLY
)
]]


cet_build_plugin(RceRawDecoder art::module
                        lardataobj_RawData
                        dune_raw_data::dune-raw-data_Overlays
			artdaq_core::artdaq-core_Data
			artdaq-core_Utilities
                        dune_raw_data::dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
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
                        BASENAME_ONLY
)

cet_build_plugin(FelixRawDecoder art::module
                        lardataobj_RawData
                        dune_raw_data::dune-raw-data_Overlays
			artdaq_core::artdaq-core_Data
			artdaq-core_Utilities
                        dune_raw_data::dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
			dune_Geometry
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
                        BASENAME_ONLY
)

cet_build_plugin(FelixIntegrityTest art::module
                        lardataobj_RawData
                        dune_raw_data::dune-raw-data_Overlays
			artdaq_core::artdaq-core_Data
			artdaq-core_Utilities
                        dune_raw_data::dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
			dune_Geometry
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
                        BASENAME_ONLY
)

cet_build_plugin(PDSPTPCRawDecoder art::module
                        lardataobj_RawData
                        dune_raw_data::dune-raw-data_Overlays
			artdaq_core::artdaq-core_Data
			artdaq-core_Utilities
                        dune_raw_data::dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
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
                        BASENAME_ONLY
)

cet_build_plugin(CRTRawDecoder art::module
                        dune_raw_data::dune-raw-data_Overlays
                        dune_Geometry
			artdaq_core::artdaq-core_Data
			artdaq-core_Utilities
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
                        BASENAME_ONLY
)

cet_build_plugin(PDSPCTBRawDecoder art::module
                        dune_raw_data::dune-raw-data_Overlays
			artdaq_core::artdaq-core_Data
			artdaq-core_Utilities
                        ART_FRAMEWORK_CORE
                        ART_FRAMEWORK_PRINCIPAL
                        ART_FRAMEWORK_SERVICES_REGISTRY
                        ART_ROOT_IO_TFILE_SUPPORT
                        ROOT::Core
                        ART_ROOT_IO_TFILESERVICE_SERVICE
                        ART_PERSISTENCY_PROVENANCE
                        messagefacility::MF_MessageLogger
                        messagefacility::MF_MessageLogger
                        CETLIB
                        ROOT_BASIC_LIB_LIST
                        RCEDAMLIB
                        BASENAME_ONLY
)

art_make_library( LIBRARY_NAME dunetpc_protodune_RawDecodingUtils
                  SOURCE
                  RawDecodingUtils.cc
                  LIBRARIES
                  lardataobj_RawData
                  dune_raw_data::dune-raw-data_Overlays
                  artdaq_core::artdaq-core_Data
		  artdaq-core_Utilities
                  ART_FRAMEWORK_CORE
                  ART_FRAMEWORK_PRINCIPAL
                  ART_FRAMEWORK_SERVICES_REGISTRY
                  ART_ROOT_IO_TFILE_SUPPORT
                  ROOT::Core
                  ART_ROOT_IO_TFILESERVICE_SERVICE
 #                 art_Persistency_Common canvas_Persistency_Common
 #                 art_Persistency_Provenance canvas_Persistency_Provenance
 #                 art_Utilities canvas_Utilities
                  messagefacility::MF_MessageLogger

                  CETLIB
                  ROOT_BASIC_LIB_LIST
                  RCEDAMLIB
)


install_headers()
install_fhicl()
install_source()
install_scripts()
