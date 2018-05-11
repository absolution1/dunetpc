# dunetpc/dune/DataPrep/Utility/CMakeLists.txt

#message(STATUS "ROOTSYS=${ROOTSYS}")
#message(STATUS "ROOT_BASIC_LIB_LIST=${ROOT_BASIC_LIB_LIST}")

art_make(BASENAME_ONLY
         LIB_LIBRARIES
         dune_ArtSupport
         dune_DuneCommon
         art_Utilities
         canvas
         cetlib
         ${ROOT_BASIC_LIB_LIST}
         NO_PLUGINS
        )

install_source()
install_headers()
