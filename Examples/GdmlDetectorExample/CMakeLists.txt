# get list of all source files
file (GLOB_RECURSE src_files "src/*.cpp" "include/*.*pp")

# define library target
add_library(ACTSGdmlDetectorExample SHARED ${src_files})
target_include_directories(ACTSGdmlDetectorExample PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_compile_options(ACTSGdmlDetectorExample PUBLIC -std=c++14)
target_link_libraries(ACTSGdmlDetectorExample PUBLIC ACTSCore ACTSTGeoPlugin)

# set installation directories
install (TARGETS ACTSGdmlDetectorExample 
	 EXPORT ACTSTargets
	 LIBRARY DESTINATION lib
	 RUNTIME DESTINATION bin)

install (DIRECTORY include/ACTS DESTINATION include)
