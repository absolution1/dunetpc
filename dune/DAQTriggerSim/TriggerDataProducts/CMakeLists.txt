# for classes_def.xml:
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )



art_make( 
          LIB_LIBRARIES lardataobj_RawData
	  art_Persistency_Common canvas_Persistency_Common
	  art_Persistency_Provenance canvas_Persistency_Provenance
	  ${MF_MESSAGELOGGER}
	  ${MF_UTILITIES}
	  ${CETLIB}
	  ${ROOT_BASIC_LIB_LIST}
	)



install_headers()
install_source()
install_fhicl()
