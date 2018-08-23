#Make sure system libraries can be found
link_directories( /usr/local/lib )

add_library(plot SHARED ChannelView.cpp TwoAxisView.cpp)
target_link_libraries(plot crtCore::persistency crtCore::CRTID)
install(TARGETS plot EXPORT crtAlgTargets DESTINATION lib)
install(FILES ChannelView.h TwoAxisView.h DESTINATION include/plot)
install(EXPORT crtAlgTargets NAMESPACE crtAlg:: DESTINATION lib/cmake/crtAlg)
