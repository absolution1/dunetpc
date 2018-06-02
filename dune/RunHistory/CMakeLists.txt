include_directories ( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories ( $ENV{NUTOOLS_DIR}/include )
include_directories( $ENV{LIBWDA_FQ_DIR}/include )
#SET(ENV{LIBWDA_LIB} "$ENV{LIBWDA_DIR}/lib")
cet_find_library( WDALIB NAMES wda PATHS ENV LIBWDA_LIB NO_DEFAULT_PATH )

FILE (GLOB src_files *.cxx)

art_make_library(  LIBRARY_NAME  RunHistoryDUNE
		             SOURCE ${src_files}
			     LIBRARIES
	                wda
	                nutools_IFDatabase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 			art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        
                        cetlib cetlib_except
			${FHICLCPP}
                        ${ROOT_BASIC_LIB_LIST}
			${Boost_SYSTEM_LIBRARY}
			lbne-raw-data_Services_ChannelMap_ChannelMapService_service
        )

art_make_exec( NAME RunInfoDUNE
               SOURCE RunInfoDUNE.cc                 
               LIBRARIES
			RunHistoryDUNE
                         #wda
			 # nutools_IFDatabase
               )

art_make_exec( NAME getDUNEPedestals
               SOURCE getDUNEPedestals.cc                 
               LIBRARIES
			RunHistoryDUNE
               )


install_headers()
install_fhicl()
install_source()

#add_subdirectory(test)
