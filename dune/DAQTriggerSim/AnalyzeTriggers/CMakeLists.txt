art_make( BASENAME_ONLY
          MODULE_LIBRARIES dune_DAQTriggerSim_Service_TemplateTriggerService_service
	                   dune_DAQTriggerSim_TriggerDataProducts
                           TriggerAlgorithms
                           larcore_Geometry_Geometry_service
                           larcorealg_Geometry
                           #lardataobj_RawData
                           nusimdata_SimulationBase	
                           ART_FRAMEWORK_CORE
                           ART_FRAMEWORK_PRINCIPAL
                           ART_FRAMEWORK_SERVICES_REGISTRY
                           ART_ROOT_IO_TFILE_SUPPORT
                           ROOT::Core
                           ART_ROOT_IO_TFILESERVICE_SERVICE
                           art_Persistency_Common canvas
                           art_Persistency_Provenance canvas
                           art_Utilities canvas
                           messagefacility::MF_MessageLogger
                           
                           cetlib cetlib_except
                           ROOT_BASIC_LIB_LIST
	 )

install_fhicl()
install_headers()
install_source()
