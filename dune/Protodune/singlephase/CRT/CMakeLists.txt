add_subdirectory(data)
add_subdirectory(alg)

include_directories( $ENV{LIBWDA_FQ_DIR}/include )
include_directories( $ENV{IFDH_ART_INC} )
include_directories( $ENV{IFDHC_FQ_DIR}/inc )
include_directories( $ENV{IFBEAM_FQ_DIR}/include )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

art_make( MODULE_LIBRARIES
                        lardataobj_RawData
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
                        lardataobj_Simulation
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
                        ${ROOT_GENVECTOR}
                        ${ROOT_GEOM}
)

install_headers()
install_fhicl()
install_source()
install_scripts()
