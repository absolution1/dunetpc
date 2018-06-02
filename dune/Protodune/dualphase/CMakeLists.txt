simple_plugin(
												Purity "module"
												dune_DetSim
                        lardataobj_RecoBase
                        larreco_Calorimetry
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_ArtDataHelper
                        lardata_Utilities
                        lardataobj_Simulation
                        larsim_MCCheater_BackTrackerService_service
           larsim_MCCheater_ParticleInventoryService_service
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
                        BASENAME_ONLY
)

simple_plugin(
			DEdxdp "module"
			dune_DetSim
                        lardataobj_RecoBase
                        larreco_Calorimetry
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_ArtDataHelper
                        lardata_Utilities
                        lardataobj_Simulation
                        larsim_MCCheater_BackTrackerService_service
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
                        BASENAME_ONLY
)

simple_plugin(
	  AnaRootParser "module"
	  larcore_Geometry_Geometry_service
	  larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
	  larsim_MCCheater_BackTrackerService_service
          larsim_MCCheater_ParticleInventoryService_service
	  lardata_Utilities
	  larevt_Filters
	  lardataobj_RawData
	  lardataobj_RecoBase
	  lardata_RecoObjects larreco_Deprecated
	  lardata_ArtDataHelper
	  larreco_RecoAlg
	  larpandora_LArPandoraInterface
	  SignalShapingServiceDUNE10kt_service
	  SignalShapingServiceDUNE35t_service
          nusimdata_SimulationBase
          ${LARRECO_LIB}
          ${LARDATA_LIB}
          ${ART_FRAMEWORK_CORE}
	  ${ART_FRAMEWORK_PRINCIPAL}
	  ${ART_FRAMEWORK_SERVICES_REGISTRY}
	  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
	  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
	  art_Persistency_Common canvas
	  art_Persistency_Provenance canvas
	  art_Utilities canvas
	  ${MF_MESSAGELOGGER}
	  
	  ${FHICLCPP}
	  cetlib cetlib_except
          ${ROOT_GEOM}
          ${ROOT_XMLIO}
          ${ROOT_GDML}
	  ${ROOT_BASIC_LIB_LIST}
)



install_headers()
install_fhicl()
install_source()
install_scripts()
