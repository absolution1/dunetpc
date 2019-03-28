simple_plugin(DetectorPropertiesServiceProtoDUNEsp "service" 
	          ProtoDUNEspDataProviders
              lardataalg_DetectorInfo
              larcore_Geometry_Geometry_service
              larcorealg_Geometry
              art_Framework_Principal
              ifdh

	      ${ART_FRAMEWORK_IO_ROOT_ROOTDB}
	      ${ART_FRAMEWORK_SERVICES_REGISTRY}
	      art_Utilities
	      canvas
	      ${SQLITE3}
	      ${MF_MESSAGELOGGER}
	      
	      cetlib cetlib_except
	      ${ROOT_BASIC_LIB_LIST}
              
)

install_headers()
install_fhicl()
install_source()