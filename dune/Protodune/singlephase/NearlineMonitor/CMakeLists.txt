
simple_plugin(SPLifetime "module"
			lardataobj_RecoBase
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
      larsim_MCCheater_BackTrackerService_service
           larsim_MCCheater_ParticleInventoryService_service
      nusimdata_SimulationBase
      ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 			art_Persistency_Common
			art_Utilities 
      canvas
			${MF_MESSAGELOGGER}
			
			cetlib cetlib_except
      ${ROOT_MATHCORE}
      ${ROOT_GENVECTOR}
      ${ROOT_GEOM}
			${ROOT_BASIC_LIB_LIST}
			BASENAME_ONLY
)

simple_plugin(PlotTrackT0 "module"
      lardataobj_RecoBase
      larcorealg_Geometry
      larcore_Geometry_Geometry_service
      larsim_MCCheater_BackTrackerService_service
           larsim_MCCheater_ParticleInventoryService_service
      nusimdata_SimulationBase
      ${ART_FRAMEWORK_CORE}
      ${ART_FRAMEWORK_PRINCIPAL}
      ${ART_FRAMEWORK_SERVICES_REGISTRY}
      ${ART_FRAMEWORK_SERVICES_OPTIONAL}
      ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
      art_Persistency_Common
      art_Utilities 
      canvas
      ${MF_MESSAGELOGGER}
      
      cetlib cetlib_except
      ${ROOT_MATHCORE}
      ${ROOT_GENVECTOR}
      ${ROOT_GEOM}
      ${ROOT_BASIC_LIB_LIST}
      BASENAME_ONLY
)

simple_plugin(PlotTrackDetails "module"
      lardataobj_RecoBase
      larcorealg_Geometry
      larcore_Geometry_Geometry_service
      larsim_MCCheater_BackTrackerService_service
           larsim_MCCheater_ParticleInventoryService_service
      nusimdata_SimulationBase
      ${ART_FRAMEWORK_CORE}
      ${ART_FRAMEWORK_PRINCIPAL}
      ${ART_FRAMEWORK_SERVICES_REGISTRY}
      ${ART_FRAMEWORK_SERVICES_OPTIONAL}
      ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
      art_Persistency_Common
      art_Utilities 
      canvas
      ${MF_MESSAGELOGGER}
      
      cetlib cetlib_except
      ${ROOT_MATHCORE}
      ${ROOT_GENVECTOR}
      ${ROOT_GEOM}
      ${ROOT_BASIC_LIB_LIST}
      BASENAME_ONLY
)

simple_plugin(PlotOpticalDetails "module"
      lardataobj_RecoBase
      ${ART_FRAMEWORK_CORE}
      ${ART_FRAMEWORK_PRINCIPAL}
      ${ART_FRAMEWORK_SERVICES_REGISTRY}
      ${ART_FRAMEWORK_SERVICES_OPTIONAL}
      ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
      art_Persistency_Common
      art_Utilities 
      canvas
      ${MF_MESSAGELOGGER}
      
      cetlib cetlib_except
      ${ROOT_MATHCORE}
      ${ROOT_GENVECTOR}
      ${ROOT_GEOM}
      ${ROOT_BASIC_LIB_LIST}
      BASENAME_ONLY
)

simple_plugin(RawEventDisplay "module"
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardataobj_RawData
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${MF_MESSAGELOGGER}
                        
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
                        BASENAME_ONLY
                        lbne_daqinput35t
)

simple_plugin(TpcMonitor "module"
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardataobj_RawData
                        nusimdata_SimulationBase
                        dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${MF_MESSAGELOGGER}
                        
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
                        BASENAME_ONLY
)

simple_plugin(SSPMonitor "module"
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardataobj_RawData
                        lardataobj_RecoBase
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${MF_MESSAGELOGGER}
                        
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
                        BASENAME_ONLY
)

simple_plugin(PDSPNearlineHeader "module"
      lardataobj_RecoBase
      ${ART_FRAMEWORK_CORE}
      ${ART_FRAMEWORK_PRINCIPAL}
      ${ART_FRAMEWORK_SERVICES_REGISTRY}
      ${ART_FRAMEWORK_SERVICES_OPTIONAL}
      ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
      art_Persistency_Common
      art_Utilities
      canvas
      ${MF_MESSAGELOGGER}
      
      cetlib cetlib_except
      ${ROOT_BASIC_LIB_LIST}
      BASENAME_ONLY
)

simple_plugin(PDSPHitMonitor "module"
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardataobj_RawData
                        nusimdata_SimulationBase
                        dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${MF_MESSAGELOGGER}

                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
                        BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
install_scripts()
