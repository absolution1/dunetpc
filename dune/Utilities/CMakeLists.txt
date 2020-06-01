#include_directories("${CMAKE_SOURCE_DIR/")

set( dune_util_lib_list lardata_Utilities_LArFFT_service
                        lardataalg_DetectorInfo
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardata_Utilities
  			${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			${ART_ROOT_IO_TFILESERVICE_SERVICE}
			${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
			art_Utilities canvas
			${MF_MESSAGELOGGER}
			
			cetlib cetlib_except
			${ROOT_BASIC_LIB_LIST}
                        ${ROOT_GEOM}
    )


simple_plugin( SignalShapingServiceDUNE  "service"
               ${dune_util_lib_list}
	       BASENAME_ONLY
        )

simple_plugin( SignalShapingServiceDUNE10kt  "service"
               ${dune_util_lib_list}
	       BASENAME_ONLY
        )

simple_plugin( SignalShapingServiceDUNEDPhase  "service"
               ${dune_util_lib_list}
	       BASENAME_ONLY
        )
	             
simple_plugin( SignalShapingServiceDUNE34kt  "service"
               ${dune_util_lib_list}
	       BASENAME_ONLY
        )
              

simple_plugin( SignalShapingServiceDUNE35t  "service"
               ${dune_util_lib_list}
	       BASENAME_ONLY
        )

simple_plugin( FileCatalogMetadataDUNE "service"
               ${dune_util_lib_list}
	       BASENAME_ONLY
	)

simple_plugin( TFileMetadataDUNE "service"
               ${dune_util_lib_list}
	       BASENAME_ONLY
	)

simple_plugin(CrpGainService   "service"
                dune_ArtSupport
                larcorealg_Geometry
		lardataobj_Simulation
		art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
		BASENAME_ONLY
             )

install_headers()
install_fhicl()
install_source()

add_subdirectory(test)
