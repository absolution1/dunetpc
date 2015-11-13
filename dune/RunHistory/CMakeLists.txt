include_directories( $ENV{LIBWDA_FQ_DIR}/include )
#SET(ENV{LIBWDA_LIB} "$ENV{LIBWDA_DIR}/lib")
cet_find_library( WDALIB NAMES wda PATHS ENV LIBWDA_LIB NO_DEFAULT_PATH )


art_make( NO_PLUGINS
          LIBRARY_NAME  RunHistoryDUNE
          LIB_LIBRARIES
	                wda
	                IFDatabase
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${CETLIB}
			${FHICLCPP}
                        ${ROOT_BASIC_LIB_LIST}
			${Boost_SYSTEM_LIBRARY}
        )



install_headers()
install_fhicl()
install_source()

#add_subdirectory(test)

