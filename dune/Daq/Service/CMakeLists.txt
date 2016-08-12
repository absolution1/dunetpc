# dunetpc/dune/Daq/Service/CMakeLists.txt
#
# David Adams
# Febraurty 2016
#
# CMake instructions for Daq services.

# Add plugin for each service.

simple_plugin(FileChannelMappingService  "service")
simple_plugin(Dune35tChannelMappingService  "service"
                lbne-raw-data_Services_ChannelMap_ChannelMapService_service
                art_Utilities canvas_Utilities
		${ROOT_CORE}
             )

add_subdirectory(test)
