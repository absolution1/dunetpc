#Make sure system libraries can be found
link_directories( /usr/local/lib )

add_library(toROOT SHARED ToFlat.cpp)
target_link_libraries(toROOT crtCore::persistency ${ROOT_LIBRARIES})
install(TARGETS toROOT EXPORT crtAlgTargets DESTINATION lib)
install(FILES ToFlat.h DESTINATION include/toROOT)
install(EXPORT crtAlgTargets NAMESPACE crtAlg:: DESTINATION lib/cmake/crtAlg)
