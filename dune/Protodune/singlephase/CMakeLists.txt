
art_make( MODULE_LIBRARIES
                        lardataobj_RawData
			lardataobj_RecoBase
			larreco_RecoAlg
			larreco_Calorimetry
			larreco_RecoAlg_PMAlg
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
                        lardata_ArtDataHelper
			lardata_Utilities
			larsim_Simulation 
                        lardataobj_Simulation
                        larsim_MCCheater_BackTrackerService_service
                        larsim_MCCheater_ParticleInventoryService_service
                        nutools_ParticleNavigation 
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
)

<<<<<<< HEAD
#simple_plugin(HadCal "module"
#                        lardataobj_RawData
#                        lardataobj_RecoBase
#                        larreco_RecoAlg
#                        larreco_Calorimetry
#                        larreco_RecoAlg_PMAlg
#                        larcorealg_Geometry
#                        larcore_Geometry_Geometry_service
#                        lardata_ArtDataHelper
#                        lardata_Utilities
#                        larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
#                        
#                        nusimdata_SimulationBase
#                        ${ART_FRAMEWORK_CORE}
#                        ${ART_FRAMEWORK_PRINCIPAL}
#                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
#                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
#                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
#                        art_Persistency_Common canvas
#                        art_Persistency_Provenance canvas
#                        art_Utilities canvas
#                        ${MF_MESSAGELOGGER}
#                        ${MF_UTILITIES}
#                        cetlib cetlib_except
#                        ${ROOT_BASIC_LIB_LIST}
#                        BASENAME_ONLY
#)
#
#simple_plugin(RecoEff "module"
#                        lardataobj_RecoBase
#                        larcorealg_Geometry
#                        larcore_Geometry_Geometry_service
#                        lardata_Utilities
#                        lardataobj_Simulation
#                        larsim_MCCheater_BackTrackerService_service
#           larsim_MCCheater_ParticleInventoryService_service
#                        ${ART_FRAMEWORK_CORE}
#                        ${ART_FRAMEWORK_PRINCIPAL}
#                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
#                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
#                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
#                        art_Persistency_Common canvas
#                        art_Persistency_Provenance canvas
#                        art_Utilities canvas
#                        ${MF_MESSAGELOGGER}
#                        ${MF_UTILITIES}
#                        cetlib cetlib_except
#                        ${ROOT_BASIC_LIB_LIST}
#                        BASENAME_ONLY
#)
#
#simple_plugin(RecoStats "module"
#                        lardataobj_RecoBase
#                        larreco_Calorimetry
#                        larcorealg_Geometry
#                        larcore_Geometry_Geometry_service
#                        lardata_ArtDataHelper
#                        lardata_Utilities
#                        lardataobj_Simulation
#                        larsim_MCCheater_BackTrackerService_service
#           larsim_MCCheater_ParticleInventoryService_service
#                        ${ART_FRAMEWORK_CORE}
#                        ${ART_FRAMEWORK_PRINCIPAL}
#                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
#                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
#                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
#                        art_Persistency_Common canvas
#                        art_Persistency_Provenance canvas
#                        art_Utilities canvas
#                        ${MF_MESSAGELOGGER}
#                        ${MF_UTILITIES}
#                        cetlib cetlib_except
#                        ${ROOT_BASIC_LIB_LIST}
#                        BASENAME_ONLY
#)
#
#simple_plugin(RecoHitsEff "module"
#                        lardataobj_RecoBase
#                        larreco_Calorimetry
#                        larcorealg_Geometry
#                        larcore_Geometry_Geometry_service
#                        lardata_ArtDataHelper
#                        lardata_Utilities
#                        lardataobj_Simulation
#                        larsim_MCCheater_BackTrackerService_service
#           larsim_MCCheater_ParticleInventoryService_service
#                        ${ART_FRAMEWORK_CORE}
#                        ${ART_FRAMEWORK_PRINCIPAL}
#                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
#                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
#                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
#                        art_Persistency_Common canvas
#                        art_Persistency_Provenance canvas
#                        art_Utilities canvas
#                        ${MF_MESSAGELOGGER}
#                        ${MF_UTILITIES}
#                        cetlib cetlib_except
#                        ${ROOT_BASIC_LIB_LIST}
#                        BASENAME_ONLY
#)
#
#simple_plugin(DEdx "module"
#                        lardataobj_RecoBase
#                        larreco_Calorimetry
#                        larcorealg_Geometry
#                        larcore_Geometry_Geometry_service
#                        lardata_ArtDataHelper
#                        lardata_Utilities
#                        lardataobj_Simulation
#                        larsim_MCCheater_BackTrackerService_service
#           larsim_MCCheater_ParticleInventoryService_service
#                        ${ART_FRAMEWORK_CORE}
#                        ${ART_FRAMEWORK_PRINCIPAL}
#                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
#                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
#                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
#                        art_Persistency_Common canvas
#                        art_Persistency_Provenance canvas
#                        art_Utilities canvas
#                        ${MF_MESSAGELOGGER}
#                        ${MF_UTILITIES}
#                        cetlib cetlib_except
#                        ${ROOT_BASIC_LIB_LIST}
#                        BASENAME_ONLY
#)
#
#simple_plugin(calcuttjRecoEff "module"
#			lardataobj_RawData
#			lardataobj_RecoBase
#			larreco_RecoAlg
#			larreco_Calorimetry
#			larreco_RecoAlg_PMAlg
#			larcorealg_Geometry
#			larcore_Geometry_Geometry_service
#			lardata_Utilities
#			larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
#                        larsim_MCCheater_BackTrackerService_service
#		        larsim_MCCheater_ParticleInventoryService_service	
#			nusimdata_SimulationBase
#			${ART_FRAMEWORK_CORE}
#			${ART_FRAMEWORK_PRINCIPAL}
#			${ART_FRAMEWORK_SERVICES_REGISTRY}
#			${ART_FRAMEWORK_SERVICES_OPTIONAL}
#			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
# 			art_Persistency_Common canvas
#			art_Persistency_Provenance canvas
#			art_Utilities canvas
#			${MF_MESSAGELOGGER}
#			${MF_UTILITIES}
#			cetlib cetlib_except
#			${ROOT_BASIC_LIB_LIST}
#			BASENAME_ONLY
#)

install_headers()
install_fhicl()
install_source()
install_scripts()
add_subdirectory(NearlineMonitor)
add_subdirectory(RawDecoding)
#add_subdirectory(calcuttjRecoEff)
