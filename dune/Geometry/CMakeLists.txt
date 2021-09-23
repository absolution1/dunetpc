#MESSAGE( STATUS "Boost_SYSTEM_LIBRARY:         " ${Boost_SYSTEM_LIBRARY} )

art_make( BASENAME_ONLY
	  EXCLUDE CheckGeometry_module.cc CheckDPhaseGeometry_module.cc issue19191_module.cc
	  LIBRARY_NAME dune_Geometry
	  LIB_LIBRARIES larcorealg_Geometry
                        messagefacility::MF_MessageLogger
                        cetlib cetlib_except

                        FHICLCPP
                        art_Framework_Services_Registry
                        art_Utilities
                        canvas
                        ROOT_BASIC_LIB_LIST
                        ROOT_GEOM
                        ${Boost_SYSTEM_LIBRARY}
          SERVICE_LIBRARIES larcorealg_Geometry
                        dune_Geometry

                        FHICLCPP
                        cetlib cetlib_except
                        ROOT_GEOM
                        ROOT_XMLIO
                        ROOT_GDML
                        ROOT_BASIC_LIB_LIST
        )

cet_build_plugin ( CheckGeometry art::module
                     larcorealg_Geometry
                     larcore_Geometry_Geometry_service
                     ART_FRAMEWORK_CORE
                     ART_FRAMEWORK_PRINCIPAL
                     ART_FRAMEWORK_SERVICES_REGISTRY
                     ART_ROOT_IO_TFILE_SUPPORT
                     ROOT::Core
                     ART_ROOT_IO_TFILESERVICE_SERVICE
                     ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE
                     art_Persistency_Common canvas
                     art_Persistency_Provenance canvas
                     art_Utilities canvas
                     messagefacility::MF_MessageLogger

                     FHICLCPP
                     cetlib cetlib_except
                     ROOT_BASIC_LIB_LIST
                     ROOT_GEOM


               BASENAME_ONLY
          )

cet_build_plugin ( CheckDPhaseGeometry art::module
                     larcorealg_Geometry
                     larcore_Geometry_Geometry_service
                     ART_FRAMEWORK_CORE
                     ART_FRAMEWORK_PRINCIPAL
                     ART_FRAMEWORK_SERVICES_REGISTRY
                     ART_ROOT_IO_TFILE_SUPPORT
                     ROOT::Core
                     ART_ROOT_IO_TFILESERVICE_SERVICE
                     ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE
                     art_Persistency_Common canvas
                     art_Persistency_Provenance canvas
                     art_Utilities canvas
                     messagefacility::MF_MessageLogger

                     FHICLCPP
                     cetlib cetlib_except
                     ROOT_BASIC_LIB_LIST
                     ROOT_GEOM

               BASENAME_ONLY
          )

add_subdirectory(gdml)
add_subdirectory(test)

install_headers()
install_fhicl()
install_source()
