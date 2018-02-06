# dunetpc/dune/DuneCommon/CMakeLists.txt

ADD_SUBDIRECTORY(test)
add_subdirectory(Tool)
ADD_SUBDIRECTORY(Service)

#message(STATUS "ROOTSYS=${ROOTSYS}")
#message(STATUS "ROOT_BASIC_LIB_LIST=${ROOT_BASIC_LIB_LIST}")

art_make(BASENAME_ONLY
        LIB_LIBRARIES
        art_Persistency_Provenance canvas
        ${ROOT_BASIC_LIB_LIST}
        ${ROOTSYS}/lib/libHistPainter.so
        NO_PLUGINS
        )

install_fhicl()
install_source()
install_headers()
