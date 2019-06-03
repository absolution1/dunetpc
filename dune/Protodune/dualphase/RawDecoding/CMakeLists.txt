simple_plugin(AmcRawDecoder "module"
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



install_headers()
install_fhicl()
install_source()
install_scripts()
add_subdirectory(data)
