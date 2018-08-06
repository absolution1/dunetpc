art_make(BASENAME_ONLY
	 MODULE_LIBRARIES
	 lardataobj_RecoBase
         lardataobj_RawData
         lardata_ArtDataHelper
	 larcore_Geometry_Geometry_service
         larcorealg_Geometry
         nusimdata_SimulationBase
	 larsim_MCCheater_BackTrackerService_service
           larsim_MCCheater_ParticleInventoryService_service
	 HitFinderDUNE
         ${ART_FRAMEWORK_CORE}
	 ${ART_FRAMEWORK_PRINCIPAL}
	 ${ART_FRAMEWORK_SERVICES_REGISTRY}
         ${ART_FRAMEWORK_SERVICES_OPTIONAL}
         ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 	 art_Persistency_Common canvas
	 art_Persistency_Provenance canvas
	 art_Utilities canvas
	 ${MF_MESSAGELOGGER}
	 
	 cetlib cetlib_except
	 ${ROOT_BASIC_LIB_LIST}
	 lbne_daqinput35t
         ${Boost_SYSTEM_LIBRARY}
)

install_headers()
install_fhicl()
install_scripts()
