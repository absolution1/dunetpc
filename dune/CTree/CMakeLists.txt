art_make( BASENAME_ONLY
          MODULE_LIBRARIES Geometry
          dune_Geometry
          Geometry_service
          Simulation
          BackTracker_service
          Utilities
          RawData
          RecoBase
          RecoObjects
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
