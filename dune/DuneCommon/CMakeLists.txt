# dunetpc/dune/DuneCommon/CMakeLists.txt

ADD_SUBDIRECTORY(test)
add_subdirectory(Tool)
ADD_SUBDIRECTORY(Service)

art_make(BASENAME_ONLY
        LIB_LIBRARIES
        art_Persistency_Provenance canvas
        ${ROOT_BASIC_LIB_LIST}
        NO_PLUGINS
        )

install_fhicl()
install_source()
install_headers()
