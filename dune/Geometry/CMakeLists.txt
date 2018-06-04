#MESSAGE( STATUS "Boost_SYSTEM_LIBRARY:         " ${Boost_SYSTEM_LIBRARY} )

art_make( BASENAME_ONLY
	  EXCLUDE CheckGeometry_module.cc CheckDPhaseGeometry_module.cc
	  LIBRARY_NAME dune_Geometry
	  LIB_LIBRARIES larcorealg_Geometry
                        ${MF_MESSAGELOGGER}
                        cetlib cetlib_except
                        ${FHICLCPP}
                        art_Framework_Services_Registry
                        art_Utilities
                        canvas
                        ${ROOT_BASIC_LIB_LIST}
                        ${ROOT_GEOM}
                        ${Boost_SYSTEM_LIBRARY}
	  SERVICE_LIBRARIES larcorealg_Geometry
	  		dune_Geometry
			
	                ${FHICLCPP}
			cetlib cetlib_except
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
			${ROOT_BASIC_LIB_LIST}
        )							   

simple_plugin ( CheckGeometry "module"
                     larcorealg_Geometry
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
                     
                     ${FHICLCPP}
                     cetlib cetlib_except
                     ${ROOT_BASIC_LIB_LIST}
                     ${ROOT_GEOM}


               BASENAME_ONLY
          )

simple_plugin ( CheckDPhaseGeometry "module"
                     larcorealg_Geometry
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
