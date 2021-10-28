add_subdirectory(data)
add_subdirectory(alg)

include_directories( $ENV{LIBWDA_FQ_DIR}/include )
include_directories( $ENV{IFDH_ART_INC} )
include_directories( $ENV{IFDHC_FQ_DIR}/inc )
include_directories( $ENV{IFBEAM_FQ_DIR}/include )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

# temporary fix for an extra semicolon in the root v6_22_06a TVirtualPaveStats header
cet_add_compiler_flags(CXX -Wno-pedantic)

art_make( MODULE_LIBRARIES
                        lardataalg_DetectorInfo
                        lardataobj_RawData
                        lardataalg_DetectorInfo
			larcorealg_Geometry
                        larcorealg_CoreUtils
			larcore_Geometry_Geometry_service
                        lardataobj_Simulation
                        nusimdata_SimulationBase
			lardataobj_RecoBase
			lardataobj_AnalysisBase
			larsim_MCCheater_ParticleInventoryService_service
			larsim_MCCheater_BackTrackerService_service
			larsim_MCCheater_PhotonBackTrackerService_service
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
                        ROOT_GENVECTOR
                        ROOT_GEOM
)

install_headers()
install_fhicl()
install_source()
install_scripts()
