cet_enable_asserts()
art_make( NO_PLUGINS
          LIBRARY_NAME  ProtoDUNEspDataProviders
          LIB_LIBRARIES
                        lardataalg_DetectorInfo
  	                wda
	                nuevdb_IFDatabase
                        larcorealg_Geometry
                        art_Framework_Core
                        art_Framework_Principal
                        ART_ROOT_IO_ROOTDB
                        art_Utilities
                        messagefacility::MF_MessageLogger
                        ART_FRAMEWORK_SERVICES_REGISTRY
                        cetlib cetlib_except
                        canvas
                        ifdh
                        FHICLCPP
                        ROOT_BASIC_LIB_LIST
                        ${Boost_SYSTEM_LIBRARY}
        )


install_headers()
install_fhicl()
install_source()
