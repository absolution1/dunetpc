# dune/ArtSupport/CMakeLists.txt
#
# David Adams
# September 2015
#
# CMake instructions for the AXServices library which holds
# helper classes for dealing with art services.

art_make(
  LIB_LIBRARIES
    art_Framework_EventProcessor
    art_Framework_Services_Registry
    art_Utilities canvas
    fhiclcpp
    cetlib cetlib_except
    ${Boost_SYSTEM_LIBRARY}
    art_Framework_Services_System_TriggerNamesService_service
    art_Framework_Services_System_CurrentModule_service
    ${ROOT_CORE}
)

add_subdirectory(Tool)
add_subdirectory(test)
add_subdirectory(exe)

install_headers()
install_fhicl()
install_source( README.txt )
