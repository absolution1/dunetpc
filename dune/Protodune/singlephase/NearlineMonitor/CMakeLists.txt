
simple_plugin(SPLifetime "module"
			lardataobj_RecoBase
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
      larsim_MCCheater_BackTracker_service
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
			${MF_UTILITIES}
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
      larsim_MCCheater_BackTracker_service
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
      ${MF_UTILITIES}
      cetlib cetlib_except
      ${ROOT_MATHCORE}
      ${ROOT_GENVECTOR}
      ${ROOT_GEOM}
      ${ROOT_BASIC_LIB_LIST}
      BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
install_scripts()
