<<<<<<< HEAD
art_make(	
	  LIB_LIBRARIES
			${ART_FRAMEWORK_CORE}
			art_Utilities canvas
			${MF_MESSAGELOGGER}
			
			cetlib cetlib_except
                        ${Boost_SYSTEM_LIBRARY}
			${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES 
            lardataobj_RawData
            larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
		    larsim_MCCheater_BackTrackerService_service
            larsim_MCCheater_ParticleInventoryService_service
			${ART_FRAMEWORK_SERVICES_BASIC}
			 
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardataobj_RecoBase
			lardataobj_AnalysisBase	
			lardata_RecoObjects larreco_Deprecated  
			larreco_RecoAlg
			lardata_Utilities
			larevt_Filters
			nusimdata_SimulationBase
			${LARRECO_LIB}
			${LARDATA_LIB}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
			${MF_MESSAGELOGGER}
			
	                ${FHICLCPP}
			cetlib cetlib_except
			${CLHEP}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
			${ROOT_BASIC_LIB_LIST}
        )
=======
art_make( BASENAME_ONLY
          MODULE_LIBRARIES lardataobj_RecoBase
                           lardataobj_AnalysisBase
                           larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           lardataobj_RawData
                           lardata_Utilities
                           larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                           larsim_MCCheater_BackTrackerService_service
                           larsim_MCCheater_ParticleInventoryService_service
                           nusimdata_SimulationBase
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_FRAMEWORK_SERVICES_BASIC}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                           ${ART_PERSISTENCY_PROVENANCE}
                           art_Persistency_Common
                           art_Utilities canvas
                           ${MF_MESSAGELOGGER}
                           ${FHICLCPP}
                           cetlib cetlib_except
                           ${ROOT_BASIC_LIB_LIST}
         )
>>>>>>> tags/v08_09_00

install_headers()
install_fhicl()
install_source()

