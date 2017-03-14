art_make( BASENAME_ONLY
	  MODULE_LIBRARIES lardataobj_RawData
	  	  	   lardataobj_RecoBase	  
	  		   larcore_Geometry  
	  		   larcore_Geometry_Geometry_service
	 		   lardata_Utilities
	  		   nusimdata_SimulationBase
          		   ${ART_FRAMEWORK_CORE}
	  		   ${ART_FRAMEWORK_PRINCIPAL}
	  		   ${ART_FRAMEWORK_SERVICES_REGISTRY}
          		   ${ART_FRAMEWORK_SERVICES_OPTIONAL}
          		   ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 	  		   art_Persistency_Common canvas_Persistency_Common
	  		   art_Persistency_Provenance canvas_Persistency_Provenance
	  		   art_Utilities canvas_Utilities
	  		   ${MF_MESSAGELOGGER}
	  		   ${MF_UTILITIES}
	  		   ${CETLIB}
	  		   ${ROOT_BASIC_LIB_LIST}
			   ${ROOT_EG}
        )

install_headers()
install_fhicl()
install_source()
install_scripts()
