make_dune_tool_builder(dune::RunningSumTPFinderTool)

art_make( BASENAME_ONLY MODULE_LIBRARIES
  dune_raw_data::dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
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
  ART_FRAMEWORK_CORE
  ART_FRAMEWORK_PRINCIPAL
  ART_FRAMEWORK_SERVICES_REGISTRY
  ART_ROOT_IO_TFILE_SUPPORT
  ROOT::Core
  ART_ROOT_IO_TFILESERVICE_SERVICE
  art_Persistency_Common canvas
  art_Persistency_Provenance canvas
  art_Utilities canvas
  messagefacility::MF_MessageLogger
  
  cetlib cetlib_except
  ROOT_BASIC_LIB_LIST
  EXCLUDE read_packed.cxx RunningSumTPFinderPass1_tool.cc RunningSumTPFinderPass2_tool.cc RunningSumTPFinderPass3_tool.cc RunningSumTPFinderPass4_tool.cc
  )

cet_build_plugin(RunningSumTPFinderPass1 art::tool
  fhiclcpp cetlib cetlib_except
  ${Boost_SYSTEM_LIBRARY}
)
cet_build_plugin(RunningSumTPFinderPass2 art::tool
  fhiclcpp cetlib cetlib_except
  ${Boost_SYSTEM_LIBRARY}
)
cet_build_plugin(RunningSumTPFinderPass3 art::tool
  fhiclcpp cetlib cetlib_except
  ${Boost_SYSTEM_LIBRARY}
)
cet_build_plugin(RunningSumTPFinderPass4 art::tool
  fhiclcpp cetlib cetlib_except
  ${Boost_SYSTEM_LIBRARY}
)

install_fhicl()
install_headers()
install_source()
