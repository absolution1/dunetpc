
simple_plugin(MichelHitCounter "module"
              lardataobj_RecoBase

              larreco_Calorimetry
              larreco_RecoAlg_PMAlg

              larsim_MCCheater_BackTracker_service

              nusimdata_SimulationBase
              ${ART_FRAMEWORK_CORE}
              ${ART_FRAMEWORK_PRINCIPAL}
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
              ${MF_MESSAGELOGGER}
              ${MF_UTILITIES}
              ${CETLIB}
              ${ROOT_BASIC_LIB_LIST}
              BASENAME_ONLY

)

simple_plugin(MichelEventSelection "module"
              larcorealg_Geometry
              lardataobj_RecoBase
              larreco_Calorimetry
              larreco_RecoAlg_PMAlg
              larsim_MCCheater_BackTracker_service
              nusimdata_SimulationBase
              ${ART_FRAMEWORK_CORE}
              ${ART_FRAMEWORK_PRINCIPAL}
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
              ${MF_MESSAGELOGGER}
              ${MF_UTILITIES}
              ${CETLIB}
              ${ROOT_BASIC_LIB_LIST}
              BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
install_scripts()
