add_subdirectory(SNAnaClustering)
add_subdirectory(RunningSumHitFinder)

art_make( BASENAME_ONLY MODULE_LIBRARIES
  dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
  larcore_Geometry_Geometry_service
  lardataalg_DetectorInfo
  larcorealg_Geometry
  lardataobj_RawData
  lardataobj_RecoBase
  lardataobj_AnalysisBase
  dune_DuneObj
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  larsim_MCCheater_PhotonBackTrackerService_service
  larsim_Simulation nug4_ParticleNavigation lardataobj_Simulation
  lardata_ArtDataHelper
  nusimdata_SimulationBase	
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_BASIC}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  art_Persistency_Common canvas
  art_Persistency_Provenance canvas
  art_Utilities canvas
  ${MF_MESSAGELOGGER}
  
  cetlib cetlib_except
  ${ROOT_BASIC_LIB_LIST}
  EXCLUDE read_packed.cxx TriggerPrimitiveFinderPass1_tool.cc TriggerPrimitiveFinderPass2_tool.cc TriggerPrimitiveFinderTemplate_tool.cc
  )

simple_plugin(TriggerPrimitiveFinderPass1 "tool"
  fhiclcpp cetlib cetlib_except
  ${Boost_SYSTEM_LIBRARY}
)

simple_plugin(TriggerPrimitiveFinderPass2 "tool"
  fhiclcpp cetlib cetlib_except
  ${Boost_SYSTEM_LIBRARY}
  dune_DAQSimAna_TriggerPrimitiveFinderPass1_tool
)

simple_plugin(TriggerPrimitiveFinderTemplate "tool"
  fhiclcpp cetlib cetlib_except
  ${Boost_SYSTEM_LIBRARY}
)

cet_make_exec(read_packed SOURCE read_packed.cxx)

install_fhicl()
install_headers()
install_source()
