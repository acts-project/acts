# get list of all source files
file (GLOB src_files "src/*.cxx")

# define library target
add_library(ACTSGenericDetectorExample SHARED ${src_files})
target_include_directories(ACTSGenericDetectorExample PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_compile_options(ACTSGenericDetectorExample PUBLIC -std=c++1y)
target_link_libraries(ACTSGenericDetectorExample PUBLIC ACTSCore)

# define executables
add_executable(ACTSGenericDetector bin/main.cxx)
target_include_directories(ACTSGenericDetector PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include/)
target_compile_options(ACTSGenericDetector PRIVATE -std=c++1y)
target_link_libraries(ACTSGenericDetector PRIVATE ACTSCore ACTSGenericDetectorExample)

# set installation directories
install (TARGETS ACTSGenericDetectorExample ACTSGenericDetector
	 EXPORT ACTSTargets
	 LIBRARY DESTINATION lib
	 RUNTIME DESTINATION bin)

install (DIRECTORY include/ACTS DESTINATION include)
