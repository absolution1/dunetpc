#include_directories("${CMAKE_SOURCE_DIR/")

set( dune_util_lib_list lardata_Utilities_LArFFT_service
                        lardataalg_DetectorInfo
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardata_Utilities
  			ART_FRAMEWORK_CORE
			ART_FRAMEWORK_PRINCIPAL
			ART_FRAMEWORK_SERVICES_REGISTRY
			art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			ART_ROOT_IO_TFILE_SUPPORT
			ROOT::Core
			ART_ROOT_IO_TFILESERVICE_SERVICE
			ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE
			art_Utilities canvas
			messagefacility::MF_MessageLogger
			
			cetlib cetlib_except
			ROOT_BASIC_LIB_LIST
                        ROOT_GEOM
    )


cet_build_plugin( SignalShapingServiceDUNE  art::service
               ${dune_util_lib_list}
	       BASENAME_ONLY
        )

cet_build_plugin( SignalShapingServiceDUNE10kt  art::service
               ${dune_util_lib_list}
	       BASENAME_ONLY
        )

cet_build_plugin( SignalShapingServiceDUNEDPhase  art::service
               ${dune_util_lib_list}
	       BASENAME_ONLY
        )
	             
cet_build_plugin( SignalShapingServiceDUNE34kt  art::service
               ${dune_util_lib_list}
	       BASENAME_ONLY
        )
              

cet_build_plugin( SignalShapingServiceDUNE35t  art::service
               ${dune_util_lib_list}
	       BASENAME_ONLY
        )

cet_build_plugin( FileCatalogMetadataDUNE art::service
               ${dune_util_lib_list}
	       BASENAME_ONLY
	)

cet_build_plugin( TFileMetadataDUNE art::service
               ${dune_util_lib_list}
	       BASENAME_ONLY
	)

cet_build_plugin(CrpGainService   art::service
                dune_ArtSupport
                larcorealg_Geometry
		lardataobj_Simulation
		art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
		BASENAME_ONLY
             )

install_headers()
install_fhicl()
install_source()

add_subdirectory(test)
