

cet_build_plugin(
           CRTGen art::module
           nusimdata_SimulationBase
           cetlib
           cetlib_except
           FHICLCPP
           CLHEP
           larcorealg_Geometry
           larcoreobj_SummaryData
           larsim_EventGenerator
           nurandom_RandomUtils_NuRandomService_service
           ART_FRAMEWORK_SERVICES_REGISTRY
           ART_ROOT_IO_TFILESERVICE_SERVICE
           ART_ROOT_IO_TFILE_SUPPORT
           CLHEP
           messagefacility::MF_MessageLogger
  	   ROOT_BASIC_LIB_LIST
           ROOT_GEOM
           ROOT_XMLIO
           ROOT_GDML
)
if(FALSE)
cet_build_plugin(
			Purity art::module
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
                        BASENAME_ONLY
)
endif()
cet_build_plugin(
			DEdxdp art::module
			dune_DetSim
                        lardataobj_RecoBase
                        larreco_Calorimetry
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_ArtDataHelper
                        lardata_Utilities
                        lardataobj_Simulation
                        larsim_MCCheater_BackTrackerService_service
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
                        BASENAME_ONLY
)

cet_build_plugin(
	  AnaRootParser art::module
	  larcore_Geometry_Geometry_service
	  larsim_Simulation nug4_ParticleNavigation lardataobj_Simulation
	  larsim_MCCheater_BackTrackerService_service
          larsim_MCCheater_ParticleInventoryService_service
	  lardata_Utilities
	  larevt_Filters
	  lardataobj_RawData
	  lardataobj_RecoBase
	  lardata_RecoObjects
	  lardata_ArtDataHelper
	  lardataalg_DetectorInfo
	  larcorealg_Geometry
	  larcorealg_GeoAlgo
	  larreco_RecoAlg
	  larpandora_LArPandoraInterface
	  SignalShapingServiceDUNE10kt_service
	  SignalShapingServiceDUNE35t_service
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
	  
	  FHICLCPP
	  cetlib cetlib_except
          ROOT_GEOM
          ROOT_XMLIO
          ROOT_GDML
	  ROOT_BASIC_LIB_LIST
)



install_headers()
install_fhicl()
install_source()
install_scripts()

add_subdirectory(Light)
add_subdirectory(RawDecoding)
add_subdirectory(fcl)
#add_subdirectory(test)
