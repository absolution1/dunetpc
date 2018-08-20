include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

art_make( LIB_LIBRARIES
                        lardataobj_Simulation
                        nusimdata_SimulationBase
			${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
 			art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
			cetlib cetlib_except
			${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_source()
