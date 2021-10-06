# dunetpc/dune/DuneCommon/CMakeLists.txt

add_subdirectory(test)

#message(STATUS "ROOTSYS=${ROOTSYS}")
#message(STATUS "ROOT_BASIC_LIB_LIST=${ROOT_BASIC_LIB_LIST}")

art_make(BASENAME_ONLY
         dune_DuneExample_dict
         LIB_LIBRARIES
         art_Persistency_Provenance canvas
         ROOT_BASIC_LIB_LIST
         ROOT::HistPainter
         NO_PLUGINS
        )

#install_fhicl()
install_source()
install_headers()
