# dunetpc/dune/DuneCommon/Utility/CMakeLists.txt

add_subdirectory(test)

#message(STATUS "ROOTSYS=${ROOTSYS}")
#message(STATUS "ROOT_BASIC_LIB_LIST=${ROOT_BASIC_LIB_LIST}")

# Code to display all cmake variables
#get_cmake_property(_variableNames VARIABLES)
#list (SORT _variableNames)
#foreach (_variableName ${_variableNames})
#  message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

art_make(BASENAME_ONLY
         LIB_LIBRARIES
         art_Persistency_Provenance
         canvas
         ${ROOT_BASIC_LIB_LIST}
         ${ROOTSYS}/lib/libHistPainter.so
         $ENV{FFTW_LIBRARY}/libfftw3.so
         NO_PLUGINS
        )

install_source()
install_headers()
