# dunetpc/dune/DuneCommon/Service/CMakeLists.txt

# David Adams
# October 2016

# Add plugin for each service.

cet_build_plugin(GeoRopChannelGroupService   art::service
                larcorealg_Geometry
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(GeoApaChannelGroupService   art::service
                larcorealg_Geometry
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(FixedChannelGroupService   art::service
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

add_subdirectory(test)
