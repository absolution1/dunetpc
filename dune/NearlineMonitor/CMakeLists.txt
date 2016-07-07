add_subdirectory(evd)
add_subdirectory(channelMap)

simple_plugin(NearlineAna "module"
                           lbne-raw-data_Overlays
                           lbne_daqinput35t       #For online to offline channel map
                           larcore_Geometry_Geometry_service
                           larcore_Geometry
                           lardata_RawData
                           ${SIMULATIONBASE}
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                           ${ART_PERSISTENCY_COMMON}
                           ${ART_PERSISTENCY_PROVENANCE}
                           ${ART_UTILITIES}
                           ${MF_MESSAGELOGGER}
                           ${MF_UTILITIES}
                           ${CETLIB}
                           ${ROOT_BASIC_LIB_LIST}
                           BASENAME_ONLY
)

art_make_exec( NearlineSummariser  
               SOURCE NearlineSummariser.cxx
               LIBRARIES ${ROOT_BASIC_LIB_LIST}
)

install_fhicl()
install_headers()

