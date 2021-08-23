cet_build_plugin(SPLifetime art::module
              lardataobj_RecoBase
              larcorealg_Geometry
              larcore_Geometry_Geometry_service
              larsim_MCCheater_BackTrackerService_service
              larsim_MCCheater_ParticleInventoryService_service
              nusimdata_SimulationBase
              ART_FRAMEWORK_CORE
              ART_FRAMEWORK_PRINCIPAL
              ART_FRAMEWORK_SERVICES_REGISTRY
              ART_ROOT_IO_TFILE_SUPPORT
              ROOT::Core
              ART_ROOT_IO_TFILESERVICE_SERVICE
              ART_PERSISTENCY_PROVENANCE
              art_Persistency_Common
              art_Utilities
              canvas
              messagefacility::MF_MessageLogger
              cetlib cetlib_except
              ROOT_MATHCORE
              ROOT_GENVECTOR
              ROOT_GEOM
              ROOT_BASIC_LIB_LIST
              BASENAME_ONLY)

cet_build_plugin(PlotTrackT0 art::module
              lardataobj_RecoBase
              larcorealg_Geometry
              larcore_Geometry_Geometry_service
              larsim_MCCheater_BackTrackerService_service
              larsim_MCCheater_ParticleInventoryService_service
              nusimdata_SimulationBase
              ART_FRAMEWORK_CORE
              ART_FRAMEWORK_PRINCIPAL
              ART_FRAMEWORK_SERVICES_REGISTRY
              ART_ROOT_IO_TFILE_SUPPORT
              ROOT::Core
              ART_ROOT_IO_TFILESERVICE_SERVICE
              ART_PERSISTENCY_PROVENANCE
              art_Persistency_Common
              art_Utilities
              canvas
              messagefacility::MF_MessageLogger
              cetlib cetlib_except
              ROOT_MATHCORE
              ROOT_GENVECTOR
              ROOT_GEOM
              ROOT_BASIC_LIB_LIST
              BASENAME_ONLY)

cet_build_plugin(PlotTrackDetails art::module
              lardataobj_RecoBase
              larcorealg_Geometry
              larcore_Geometry_Geometry_service
              larsim_MCCheater_BackTrackerService_service
              larsim_MCCheater_ParticleInventoryService_service
              nusimdata_SimulationBase
              ART_FRAMEWORK_CORE
              ART_FRAMEWORK_PRINCIPAL
              ART_FRAMEWORK_SERVICES_REGISTRY
              ART_ROOT_IO_TFILE_SUPPORT
              ROOT::Core
              ART_ROOT_IO_TFILESERVICE_SERVICE
              ART_PERSISTENCY_PROVENANCE
              art_Persistency_Common
              art_Utilities
              canvas
              messagefacility::MF_MessageLogger
              cetlib cetlib_except
              ROOT_MATHCORE
              ROOT_GENVECTOR
              ROOT_GEOM
              ROOT_BASIC_LIB_LIST
              BASENAME_ONLY)

cet_build_plugin(PlotEventDetails art::module
      lardataobj_RecoBase
      larcorealg_Geometry
      larcore_Geometry_Geometry_service
      larsim_MCCheater_BackTrackerService_service
           larsim_MCCheater_ParticleInventoryService_service
      nusimdata_SimulationBase
      ART_FRAMEWORK_CORE
      ART_FRAMEWORK_PRINCIPAL
      ART_FRAMEWORK_SERVICES_REGISTRY
      ART_ROOT_IO_TFILE_SUPPORT
      ROOT::Core
      ART_ROOT_IO_TFILESERVICE_SERVICE
      ART_PERSISTENCY_PROVENANCE
      art_Persistency_Common
      art_Utilities 
      canvas
      messagefacility::MF_MessageLogger
      cetlib cetlib_except
      ROOT_MATHCORE
      ROOT_GENVECTOR
      ROOT_GEOM
      ROOT_BASIC_LIB_LIST
      BASENAME_ONLY
)

cet_build_plugin(PlotOpticalDetails art::module
              lardataobj_RecoBase
              ART_FRAMEWORK_CORE
              ART_FRAMEWORK_PRINCIPAL
              ART_FRAMEWORK_SERVICES_REGISTRY
              ART_ROOT_IO_TFILE_SUPPORT
              ROOT::Core
              ART_ROOT_IO_TFILESERVICE_SERVICE
              ART_PERSISTENCY_PROVENANCE
              art_Persistency_Common
              art_Utilities
              canvas
              messagefacility::MF_MessageLogger
              cetlib cetlib_except
              ROOT_MATHCORE
              ROOT_GENVECTOR
              ROOT_GEOM
              ROOT_BASIC_LIB_LIST
              BASENAME_ONLY)

cet_build_plugin(RawEventDisplay art::module
              larcorealg_Geometry
              larcore_Geometry_Geometry_service
              lardataalg_DetectorInfo
              lardataobj_RawData
              nusimdata_SimulationBase
              ART_FRAMEWORK_CORE
              ART_FRAMEWORK_PRINCIPAL
              ART_FRAMEWORK_SERVICES_REGISTRY
              ART_ROOT_IO_TFILE_SUPPORT
              ROOT::Core
              ART_ROOT_IO_TFILESERVICE_SERVICE
              messagefacility::MF_MessageLogger
              CETLIB
              ROOT_BASIC_LIB_LIST
              BASENAME_ONLY)

cet_build_plugin(TpcMonitor art::module
              larcorealg_Geometry
              larcore_Geometry_Geometry_service
              lardataalg_DetectorInfo
              lardataobj_RawData
              nusimdata_SimulationBase
              dune_raw_data::dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
              ART_FRAMEWORK_CORE
              ART_FRAMEWORK_PRINCIPAL
              ART_FRAMEWORK_SERVICES_REGISTRY
              ART_ROOT_IO_TFILE_SUPPORT
              ROOT::Core
              ART_ROOT_IO_TFILESERVICE_SERVICE
              ART_PERSISTENCY_PROVENANCE
              messagefacility::MF_MessageLogger
              CETLIB
              ROOT_BASIC_LIB_LIST
              BASENAME_ONLY)

cet_build_plugin(SSPMonitor art::module
              larcorealg_Geometry
              larcore_Geometry_Geometry_service
              lardataobj_RawData
              lardataobj_RecoBase
              nusimdata_SimulationBase
              ART_FRAMEWORK_CORE
              ART_FRAMEWORK_PRINCIPAL
              ART_FRAMEWORK_SERVICES_REGISTRY
              ART_ROOT_IO_TFILE_SUPPORT
              ROOT::Core
              ART_ROOT_IO_TFILESERVICE_SERVICE
              ART_PERSISTENCY_PROVENANCE
              messagefacility::MF_MessageLogger
              CETLIB
              ROOT_BASIC_LIB_LIST
              BASENAME_ONLY)

cet_build_plugin(PDSPNearlineHeader art::module
              lardataobj_RecoBase
              ART_FRAMEWORK_CORE
              ART_FRAMEWORK_PRINCIPAL
              ART_FRAMEWORK_SERVICES_REGISTRY
              ART_ROOT_IO_TFILE_SUPPORT
              ROOT::Core
              ART_ROOT_IO_TFILESERVICE_SERVICE
              art_Persistency_Common
              art_Utilities
              canvas
              messagefacility::MF_MessageLogger
              cetlib cetlib_except
              ROOT_BASIC_LIB_LIST
              BASENAME_ONLY)

cet_build_plugin(PDSPHitMonitor art::module
              larcorealg_Geometry
              larcore_Geometry_Geometry_service
              lardataobj_RawData
              nusimdata_SimulationBase
              dune_raw_data::dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
              ART_FRAMEWORK_CORE
              ART_FRAMEWORK_PRINCIPAL
              ART_FRAMEWORK_SERVICES_REGISTRY
              ART_ROOT_IO_TFILE_SUPPORT
              ROOT::Core
              ART_ROOT_IO_TFILESERVICE_SERVICE
              ART_PERSISTENCY_PROVENANCE
              messagefacility::MF_MessageLogger
              CETLIB
              ROOT_BASIC_LIB_LIST
              BASENAME_ONLY)

install_headers()
install_fhicl()
install_source()
install_scripts()
