art_make( BASENAME_ONLY
          MODULE_LIBRARIES larcore_Geometry
          dune_Geometry
          larcore_Geometry_Geometry_service
          larsim_Simulation
          larsim_MCCheater_BackTracker_service
          lardata_Utilities
          lardata_RawData
          lardata_RecoBase
          lardata_RecoObjects
          ${SIMULATIONBASE}
          ${LARRECO_LIB}
          ${LARDATA_LIB}
          ${ART_FRAMEWORK_CORE}
          ${ART_FRAMEWORK_PRINCIPAL}
          ${ART_FRAMEWORK_SERVICES_REGISTRY}
          ${ART_FRAMEWORK_SERVICES_OPTIONAL}
          ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
          ${ART_PERSISTENCY_COMMON}
          ${ART_PERSISTENCY_PROVENANCE}
          ${ART_UTILITIES}
          ${MF_MESSAGELOGGER}
          ${MF_UTILITIES}
          ${FHICLCPP}
          ${CETLIB}
          ${ROOT_BASIC_LIB_LIST}
        )

# install_headers()
install_fhicl()
install_source()
install_scripts()
