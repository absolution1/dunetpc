# dunetpc/dune/DuneCommon/CMakeLists.txt

ADD_SUBDIRECTORY(test)

art_make(BASENAME_ONLY
        LIB_LIBRARIES
        art_Persistency_Provenance canvas_Persistency_Provenance
        NO_PLUGINS
        )

install_fhicl()
install_source()
install_headers()
