# dunetpc/dune/DuneCommon/CMakeLists.txt

add_subdirectory(Utility)
add_subdirectory(Tool)
add_subdirectory(Service)
add_subdirectory(Exe)
add_subdirectory(fcl)

install(FILES valgrind_dune.sup DESTINATION ${dunetpc_inc_dir})
