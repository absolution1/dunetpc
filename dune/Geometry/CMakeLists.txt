
art_make( BASENAME_ONLY
	  EXCLUDE CheckGeometry_module.cc CheckDPhaseGeometry_module.cc
	  LIBRARY_NAME dune_Geometry
	  LIB_LIBRARIES larcore_Geometry
                        ${MF_MESSAGELOGGER}
                        cetlib cetlib_except
                        ${FHICLCPP}
                        ${ROOT_BASIC_LIB_LIST}
                        ${ROOT_GEOM}
	  SERVICE_LIBRARIES larcore_Geometry
	  		dune_Geometry
			${MF_UTILITIES}
	                ${FHICLCPP}
			cetlib cetlib_except
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
		     art_Persistency_Common canvas
		     art_Persistency_Provenance canvas
		     art_Utilities canvas
                     ${MF_MESSAGELOGGER}
                     ${MF_UTILITIES}
                     ${FHICLCPP}
                     cetlib cetlib_except
                     ${ROOT_BASIC_LIB_LIST}
                     ${ROOT_GEOM}


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
		     art_Persistency_Common canvas
		     art_Persistency_Provenance canvas
		     art_Utilities canvas
                     ${MF_MESSAGELOGGER}
                     ${MF_UTILITIES}
                     ${FHICLCPP}
                     cetlib cetlib_except
                     ${ROOT_BASIC_LIB_LIST}
                     ${ROOT_GEOM}

               BASENAME_ONLY
          )

add_subdirectory(gdml)
add_subdirectory(test)

install_headers()
install_fhicl()
install_source()
