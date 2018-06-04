# for classes_def.xml:
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

art_make( 
          LIB_LIBRARIES lardataobj_RawData
	  art_Persistency_Common canvas
	  art_Persistency_Provenance canvas
	  ${MF_MESSAGELOGGER}
	  
	  cetlib cetlib_except
	  ${ROOT_BASIC_LIB_LIST}
	)

install_headers()
install_source()
install_fhicl()
