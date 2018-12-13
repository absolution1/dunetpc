
include_directories( "$ENV{DUNEPDSPRCE_INC}" ) 
cet_find_library( RCEDAMLIB NAMES protodune-dam PATHS ENV DUNEPDSPRCE_LIB NO_DEFAULT_PATH )

simple_plugin(IcebergTPCRawDecoder "module"
                        lardataobj_RawData
                        dune-raw-data_Overlays
			artdaq-core_Data
                        dune-raw-data_Services_ChannelMap_IcebergChannelMapService_service
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 #                       art_Persistency_Common canvas_Persistency_Common
 #                       art_Persistency_Provenance canvas_Persistency_Provenance
 #                       art_Utilities canvas_Utilities
                        ${MF_MESSAGELOGGER}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
                        ${RCEDAMLIB}
			z
                        BASENAME_ONLY
)


install_headers()
install_fhicl()
install_source()
install_scripts()
