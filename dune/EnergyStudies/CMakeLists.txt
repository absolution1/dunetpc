
simple_plugin(dunezsanalysis "module"
  larcore_Geometry_Geometry_service
  larcorealg_Geometry
  lardataobj_RawData
  nusimdata_SimulationBase
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
  art_Persistency_Common canvas
  art_Persistency_Provenance canvas
  art_Utilities canvas
  ${MF_MESSAGELOGGER}
  
  cetlib cetlib_except
  ${ROOT_BASIC_LIB_LIST}
  BASENAME_ONLY
  )

simple_plugin(dEdx "module"
  larcore_Geometry_Geometry_service
  larcorealg_Geometry
  lardataobj_RawData
  lardataobj_RecoBase
  lardata_Utilities
  larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
  nusimdata_SimulationBase
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
  art_Persistency_Common canvas
  art_Persistency_Provenance canvas
  art_Utilities canvas
  ${MF_MESSAGELOGGER}
  
  cetlib cetlib_except
  ${ROOT_BASIC_LIB_LIST}
  BASENAME_ONLY
  )

simple_plugin(EMEnergyCalib "module"
  lardataobj_RawData
  lardataobj_RecoBase
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
  lardata_Utilities
  nusimdata_SimulationBase
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
  art_Persistency_Common canvas
  art_Persistency_Provenance canvas
  art_Utilities canvas
  ${MF_MESSAGELOGGER}
  
  cetlib cetlib_except
  ${ROOT_BASIC_LIB_LIST}
  BASENAME_ONLY
  )

simple_plugin(EMPi0Energy "module"
  lardataobj_RawData
  lardataobj_RecoBase
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
  lardata_Utilities
  nusimdata_SimulationBase
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
  art_Persistency_Common canvas
  art_Persistency_Provenance canvas
  art_Utilities canvas
  ${MF_MESSAGELOGGER}
  
  cetlib cetlib_except
  ${ROOT_BASIC_LIB_LIST}
  BASENAME_ONLY
  )

install_headers()
install_fhicl()
install_source()
install_scripts()
