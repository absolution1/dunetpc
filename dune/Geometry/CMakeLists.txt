
art_make( BASENAME_ONLY
	  EXCLUDE CheckGeometry_module.cc CheckDPhaseGeometry_module.cc
	  LIBRARY_NAME dune_Geometry
	  LIB_LIBRARIES larcore_Geometry
                        ${MF_MESSAGELOGGER}
                        ${CETLIB}
                        ${FHICLCPP}
                        ${ROOT_BASIC_LIB_LIST}
	  SERVICE_LIBRARIES larcore_Geometry
	  		dune_Geometry
			${MF_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
			${ROOT_BASIC_LIB_LIST}
        )							   

simple_plugin ( CheckGeometry "module"
                     larcore_Geometry
                     larcore_Geometry_Geometry_service
		     ${ART_FRAMEWORK_CORE}
   		     ${ART_FRAMEWORK_PRINCIPAL}
		     ${ART_FRAMEWORK_SERVICES_REGISTRY}
		     ${ART_FRAMEWORK_SERVICES_BASIC}
		     ${ART_FRAMEWORK_SERVICES_OPTIONAL}
		     ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
             	     ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
		     art_Persistency_Common canvas_Persistency_Common
		     art_Persistency_Provenance canvas_Persistency_Provenance
		     art_Utilities canvas_Utilities
                     ${MF_MESSAGELOGGER}
                     ${MF_UTILITIES}
                     ${FHICLCPP}
                     ${CETLIB}
                     ${ROOT_BASIC_LIB_LIST}

               BASENAME_ONLY
          )

simple_plugin ( CheckDPhaseGeometry "module"
                     larcore_Geometry
                     larcore_Geometry_Geometry_service
		     ${ART_FRAMEWORK_CORE}
   		     ${ART_FRAMEWORK_PRINCIPAL}
		     ${ART_FRAMEWORK_SERVICES_REGISTRY}
		     ${ART_FRAMEWORK_SERVICES_BASIC}
		     ${ART_FRAMEWORK_SERVICES_OPTIONAL}
		     ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
             	     ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
		     art_Persistency_Common canvas_Persistency_Common
		     art_Persistency_Provenance canvas_Persistency_Provenance
		     art_Utilities canvas_Utilities
                     ${MF_MESSAGELOGGER}
                     ${MF_UTILITIES}
                     ${FHICLCPP}
                     ${CETLIB}
                     ${ROOT_BASIC_LIB_LIST}

               BASENAME_ONLY
          )

add_subdirectory(gdml)

install_headers()
install_fhicl()
install_source()
