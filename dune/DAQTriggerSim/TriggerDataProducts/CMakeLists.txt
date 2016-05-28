# for classes_def.xml:
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )



art_make( 
          LIB_LIBRARIES lardata_RawData
	  ${ART_PERSISTENCY_COMMON}
	  ${ART_PERSISTENCY_PROVENANCE}
	  ${MF_MESSAGELOGGER}
	  ${MF_UTILITIES}
	  ${CETLIB}
	  ${ROOT_BASIC_LIB_LIST}
	)



install_headers()
install_source()
install_fhicl()
