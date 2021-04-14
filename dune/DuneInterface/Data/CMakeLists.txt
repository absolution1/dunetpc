# dunetpc/dune/DuneInterface/Data/CMakeLists.txt
#
# David Adams
# September 2017
#
# CMake instructions for DUNE service interface classes.

add_subdirectory(test)

art_make(BASENAME_ONLY
         NO_PLUGINS
        )

# Needed to build dictionaries
#art_make()

install_headers()
install_source( README.txt )
