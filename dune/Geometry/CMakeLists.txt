
art_make( BASENAME_ONLY
	  EXCLUDE CheckGeometry_module.cc CheckDPhaseGeometry_module.cc
	  LIBRARY_NAME dune_Geometry
	  LIB_LIBRARIES Geometry
                        ${MF_MESSAGELOGGER}
                        ${CETLIB}
                        ${FHICLCPP}
                        ${ROOT_BASIC_LIB_LIST}
	  SERVICE_LIBRARIES Geometry
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
                     Geometry
                     Geometry_service
		     ${ART_FRAMEWORK_CORE}
   		     ${ART_FRAMEWORK_PRINCIPAL}
		     ${ART_FRAMEWORK_SERVICES_REGISTRY}
		     ${ART_FRAMEWORK_SERVICES_BASIC}
		     ${ART_FRAMEWORK_SERVICES_OPTIONAL}
		     ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
             	     ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
		     ${ART_PERSISTENCY_COMMON}
		     ${ART_PERSISTENCY_PROVENANCE}
		     ${ART_UTILITIES}
                     ${MF_MESSAGELOGGER}
                     ${MF_UTILITIES}
                     ${FHICLCPP}
                     ${CETLIB}
                     ${ROOT_BASIC_LIB_LIST}

               BASENAME_ONLY
          )

simple_plugin ( CheckDPhaseGeometry "module"
                     Geometry
                     Geometry_service
		     ${ART_FRAMEWORK_CORE}
   		     ${ART_FRAMEWORK_PRINCIPAL}
		     ${ART_FRAMEWORK_SERVICES_REGISTRY}
		     ${ART_FRAMEWORK_SERVICES_BASIC}
		     ${ART_FRAMEWORK_SERVICES_OPTIONAL}
		     ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
             	     ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
		     ${ART_PERSISTENCY_COMMON}
		     ${ART_PERSISTENCY_PROVENANCE}
		     ${ART_UTILITIES}
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
