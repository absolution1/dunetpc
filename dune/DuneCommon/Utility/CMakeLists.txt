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
cet_find_library(FFTW NAMES fftw3 PATHS ENV FFTW_LIBRARY)
art_make(BASENAME_ONLY
         LIB_LIBRARIES
         art_Persistency_Provenance
         canvas
         ROOT::HistPainter
         ROOT_BASIC_LIB_LIST
         FFTW
         PUBLIC ROOT::Core
         NO_PLUGINS
        )

install_source()
install_headers()
