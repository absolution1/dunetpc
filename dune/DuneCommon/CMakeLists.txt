# dunetpc/dune/DuneCommon/CMakeLists.txt

ADD_SUBDIRECTORY(test)

art_make(BASENAME_ONLY
        LIB_LIBRARIES
        ${ART_PERSISTENCY_PROVENANCE}
        NO_PLUGINS
        )

install_fhicl()
install_source()
install_headers()
