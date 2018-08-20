# dunetpc/dune/DataPrep/Utility/CMakeLists.txt

#message(STATUS "ROOTSYS=${ROOTSYS}")
#message(STATUS "ROOT_BASIC_LIB_LIST=${ROOT_BASIC_LIB_LIST}")
#message(STATUS "CMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}")
#message(STATUS "HERE=${CMAKE_INSTALL_PREFIX}/${${product}_lib_dir}")
#message(STATUS "CMAKE_LIBRARY_OUTPUT_DIRECTORY=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")

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

# Install the Root dictionary in the build and product lib directories.
#install(FILES TickModTreeData_Dict_rdict.pcm DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#install(FILES TickModTreeData_Dict_rdict.pcm DESTINATION ${CMAKE_INSTALL_PREFIX}/${${product}_lib_dir})

install_source()
install_headers()

add_subdirectory(test)
