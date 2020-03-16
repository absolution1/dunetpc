# dunetpc/dune/DuneInterface/ToolCMakeLists.txt
#
# David Adams
# August 2017
#
# CMake instructions for DUNE service interface classes.

add_subdirectory(test)

art_make()

install_headers()
install_source( README.txt )
