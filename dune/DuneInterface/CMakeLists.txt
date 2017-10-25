# dunetpc/dune/DuneInterface/CMakeLists.txt
#
# David Adams
# December 2015
#
# CMake instructions for DUNE service interface classes.

add_subdirectory(Data)
add_subdirectory(Tool)

install_headers()
install_source( README.txt )
