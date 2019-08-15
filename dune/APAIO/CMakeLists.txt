art_make( BASENAME_ONLY
          MODULE_LIBRARIES lardataobj_RawData
                           larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           lardata_Utilities
                           larsim_Simulation
			   lardataobj_Simulation
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_FRAMEWORK_SERVICES_BASIC}
                           ${ART_ROOT_IO_TFILE_SUPPORT}
			   ${ROOT_CORE}
                           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                           ${ART_PERSISTENCY_PROVENANCE}
                           art_Persistency_Common
                           art_Utilities canvas
                           ${MF_MESSAGELOGGER}
                           ${FHICLCPP}
                           cetlib cetlib_except
                           ${ROOT_BASIC_LIB_LIST}
         )

install_headers()
install_fhicl()
install_source()

