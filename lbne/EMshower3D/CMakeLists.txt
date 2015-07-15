
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

simple_plugin(EMshower2D "module"
			RawData
			RecoAlg
			Shower2DAlg
			BackTracker_service
      RecoBase
			Simulation
			Geometry
	    Geometry_service
      Utilities    
			DetectorProperties_service
			LArProperties_service  
			  ${SIMULATIONBASE}
        ${ART_FRAMEWORK_CORE}
			  ${ART_FRAMEWORK_PRINCIPAL}
			  ${ART_FRAMEWORK_SERVICES_REGISTRY}
				${ART_FRAMEWORK_SERVICES_BASIC}
			  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			  ${ART_PERSISTENCY_COMMON}
			  ${ART_PERSISTENCY_PROVENANCE}
			  ${ART_UTILITIES}
			  ${MF_MESSAGELOGGER}
			  ${MF_UTILITIES}
	      ${FHICLCPP}
			  ${CETLIB}
        ${ROOT_GEOM}
        ${ROOT_XMLIO}
        ${ROOT_GDML}
			  ${ROOT_BASIC_LIB_LIST}
              BASENAME_ONLY
)

add_subdirectory(Shower2DAlg)

install_headers()
install_fhicl()
install_source()
