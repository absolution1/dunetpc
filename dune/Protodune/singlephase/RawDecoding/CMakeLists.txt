include_directories( "$ENV{DUNEPDSPRCE_INC}" )
cet_find_library( RCEDAMLIB NAMES protodune-dam PATHS ENV DUNEPDSPRCE_LIB NO_DEFAULT_PATH )

simple_plugin(PDSPTPCDataInterface   "tool"
				     art_utilities
				     canvas
				     cetlib
				     cetlib_except
                        lardataobj_RawData
                        dune-raw-data_Overlays
                        artdaq-core_Data
                        dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        ${ART_PERSISTENCY_PROVENANCE}
                        ${MF_MESSAGELOGGER}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
                        ${RCEDAMLIB}
                        z
             )


simple_plugin(TimingRawDecoder "module"
                        lardataobj_RawData
                        dune-raw-data_Overlays
                        dune_DuneObj
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        ${ART_PERSISTENCY_PROVENANCE}
                        ${MF_MESSAGELOGGER}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
                        BASENAME_ONLY
)

simple_plugin(SSPRawDecoder "module"
                        lardataobj_RawData
                        lardataobj_RecoBase
                        dune-raw-data_Overlays
                        dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
                        artdaq-core_Data
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        ${ART_PERSISTENCY_PROVENANCE}
                        ${MF_MESSAGELOGGER}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
                        BASENAME_ONLY
)

#[[
simple_plugin(PDWaveform "module"
                         lardataobj_RawData
                         lardataobj_RecoBase
                         dune-raw-data_Overlays
                         dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
                         artdaq-core_Data
                         ${ART_FRAMEWORK_CORE}
                         ${ART_FRAMEWORK_PRINCIPAL}
                         ${ART_FRAMEWORK_SERVICES_REGISTRY}
                         ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                         ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                         ${MF_MESSAGELOGGER}
			 dune_Geometry
			 ${CETLIB}
                         ${ROOT_BASIC_LIB_LIST}
			 ${ROOT_SPECTRUM}
                         BASENAME_ONLY
)
]]


simple_plugin(RceRawDecoder "module"
                        lardataobj_RawData
                        dune-raw-data_Overlays
                        artdaq-core_Data
                        dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        ${ART_PERSISTENCY_PROVENANCE}
                        ${MF_MESSAGELOGGER}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
                        ${RCEDAMLIB}
                        BASENAME_ONLY
)

simple_plugin(FelixRawDecoder "module"
                        lardataobj_RawData
                        dune-raw-data_Overlays
                        artdaq-core_Data
                        dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
			dune_Geometry
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        ${ART_PERSISTENCY_PROVENANCE}
                        ${MF_MESSAGELOGGER}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
                        ${RCEDAMLIB}
                        z
                        BASENAME_ONLY
)

simple_plugin(PDSPTPCRawDecoder "module"
                        lardataobj_RawData
                        dune-raw-data_Overlays
                        artdaq-core_Data
                        dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        ${ART_PERSISTENCY_PROVENANCE}
                        ${MF_MESSAGELOGGER}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
                        ${RCEDAMLIB}
                        z
                        BASENAME_ONLY
)

simple_plugin(CRTRawDecoder "module"
                        dune-raw-data_Overlays
                        dune_Geometry
			artdaq-core_Data
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        ${ART_PERSISTENCY_PROVENANCE}
                        ${MF_MESSAGELOGGER}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
                        ${RCEDAMLIB}
                        BASENAME_ONLY
)

simple_plugin(PDSPCTBRawDecoder "module"
                        dune-raw-data_Overlays
                        artdaq-core_Data
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        ${ART_PERSISTENCY_PROVENANCE}
                        ${MF_MESSAGELOGGER}
                        ${MF_MESSAGELOGGER}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
                        ${RCEDAMLIB}
                        BASENAME_ONLY
)

art_make_library( LIBRARY_NAME dunetpc_protodune_RawDecodingUtils
                  SOURCE
                  RawDecodingUtils.cc
                  LIBRARIES
                  lardataobj_RawData
                  dune-raw-data_Overlays
                  artdaq-core_Data
                  ${ART_FRAMEWORK_CORE}
                  ${ART_FRAMEWORK_PRINCIPAL}
                  ${ART_FRAMEWORK_SERVICES_REGISTRY}
                  ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
 #                 art_Persistency_Common canvas_Persistency_Common
 #                 art_Persistency_Provenance canvas_Persistency_Provenance
 #                 art_Utilities canvas_Utilities
                  ${MF_MESSAGELOGGER}

                  ${CETLIB}
                  ${ROOT_BASIC_LIB_LIST}
                  ${RCEDAMLIB}
)


install_headers()
install_fhicl()
install_source()
install_scripts()
add_subdirectory(data)
