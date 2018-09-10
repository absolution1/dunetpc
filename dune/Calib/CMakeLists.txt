include_directories ( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories ( $ENV{NUTOOLS_DIR}/include )
include_directories( $ENV{LIBWDA_FQ_DIR}/include )
#SET(ENV{LIBWDA_LIB} "$ENV{LIBWDA_DIR}/lib")
cet_find_library( WDALIB NAMES wda PATHS ENV LIBWDA_LIB NO_DEFAULT_PATH )

FILE (GLOB src_files *.cxx)

art_make_library(  LIBRARY_NAME  dune_Calib
		        SOURCE ${src_files}
			LIBRARIES
	                wda
	                nutools_IFDatabase
                        cetlib cetlib_except
                        ${MF_MESSAGELOGGER}
			${FHICLCPP}
                        ${ROOT_BASIC_LIB_LIST}
			${Boost_SYSTEM_LIBRARY}
        )
 
art_make_exec( NAME getXYZCalibProtoDUNE
               SOURCE getXYZCalibProtoDUNE.cc                 
               LIBRARIES
			dune_Calib
               )

art_make_exec( NAME getLinConstsProtoDUNE
               SOURCE getLinConstsProtoDUNE.cc                 
               LIBRARIES
			dune_Calib
               )

install_headers()
install_fhicl()
install_source()
