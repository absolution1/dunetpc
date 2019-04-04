include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

art_make( LIB_LIBRARIES
			${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
 			art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
			cetlib cetlib_except
			${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_source()
