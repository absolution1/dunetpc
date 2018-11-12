art_make(
          LIB_LIBRARIES
  	                wda
	                nutools_IFDatabase
                        larcorealg_Geometry
                        ${MF_MESSAGELOGGER}
                        
                        cetlib cetlib_except
                        ${FHICLCPP}
                        ${ROOT_BASIC_LIB_LIST}
                        ${Boost_SYSTEM_LIBRARY}
                      
               			${SQLITE3}
               			${ART_FRAMEWORK_IO_ROOT_ROOTDB}
)
install_fhicl()
install_source()
install_scripts()
