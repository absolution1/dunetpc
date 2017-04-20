simple_plugin(TestHardwareMapper "module"
                           dune_DAQTriggerSim_Service_HardwareMapperService_service
                     	   lbne-raw-data_Overlays
                           lbne_daqinput35t       #For online to offline channel map
                           larcore_Geometry_Geometry_service
                           larcore_Geometry
                           lardataobj_RawData
                           nusimdata_SimulationBase	
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                           art_Persistency_Common canvas
                           art_Persistency_Provenance canvas
                           art_Utilities canvas
                           ${MF_MESSAGELOGGER}
                           ${MF_UTILITIES}
                           cetlib cetlib_except
                           ${ROOT_BASIC_LIB_LIST}
                           BASENAME_ONLY
)

install_fhicl()
install_headers()
install_source()
