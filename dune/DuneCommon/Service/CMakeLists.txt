# dunetpc/dune/DuneCommon/Service/CMakeLists.txt

# David Adams
# October 2016

# Add plugin for each service.

simple_plugin(GeoRopChannelGroupService   "service"
                larcorealg_Geometry
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(GeoApaChannelGroupService   "service"
                larcorealg_Geometry
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(Dune35tChannelGroupService   "service"
                lbne-raw-data_Services_ChannelMap_ChannelMapService_service
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(FixedChannelGroupService   "service"
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

add_subdirectory(test)

