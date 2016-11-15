# dunetpc/dune/DuneCommon/Service/CMakeLists.txt

# David Adams
# October 2016

# Add plugin for each service.

simple_plugin(GeoRopChannelGroupService   "service"
                larcore_Geometry
                art_Utilities canvas_Utilities
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(GeoApaChannelGroupService   "service"
                larcore_Geometry
                art_Utilities canvas_Utilities
                ${ROOT_BASIC_LIB_LIST}
             )

add_subdirectory(test)

