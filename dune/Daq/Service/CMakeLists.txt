# dunetpc/dune/Daq/Service/CMakeLists.txt
#
# David Adams
# Febraurty 2016
#
# CMake instructions for Daq services.

# Add plugin for each service.

cet_build_plugin(FileChannelMappingService  art::service cetlib_except)

add_subdirectory(test)
