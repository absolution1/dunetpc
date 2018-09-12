include_directories( $ENV{LIBWDA_FQ_DIR}/include )
include_directories( $ENV{IFDH_ART_INC} )
include_directories( $ENV{IFDHC_FQ_DIR}/inc )
include_directories( $ENV{IFBEAM_FQ_DIR}/include )

cet_find_library( WDALIB NAMES wda PATHS ENV LIBWDA_LIB NO_DEFAULT_PATH )
cet_find_library( IFBEAMSERVICE NAMES IFBeam_service PATHS ENV IFDH_ART_LIB NO_DEFAULT_PATH )
cet_find_library( IFBEAM NAMES ifbeam PATHS ENV IFBEAM_LIB   NO_DEFAULT_PATH )


art_make( MODULE_LIBRARIES
			${IFBEAM}
			${IFBEAMSERVICE}
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
                        dune_DuneObj
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
			
			cetlib cetlib_except
			${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_fhicl()
install_source()
install_scripts()
add_subdirectory(NearlineMonitor)
add_subdirectory(RawDecoding)
#add_subdirectory(calcuttjRecoEff)
add_subdirectory(CRT)
add_subdirectory(CTB)
add_subdirectory(Tool)
add_subdirectory(dEdxcalibration)
add_subdirectory(XYZcalibration)
