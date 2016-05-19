# get list of all source files
file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp")

# define library target
add_library(ACTSGenericDetectorExample SHARED ${src_files})
target_include_directories(ACTSGenericDetectorExample PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_compile_options(ACTSGenericDetectorExample PUBLIC -std=c++14)
target_link_libraries(ACTSGenericDetectorExample PUBLIC ACTSCore)

# define executables
add_executable(ACTSGenericDetector bin/main.cxx)
target_include_directories(ACTSGenericDetector PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include/)
target_compile_options(ACTSGenericDetector PRIVATE -std=c++14)
target_link_libraries(ACTSGenericDetector PRIVATE ACTSCore ACTSGenericDetectorExample)

# set installation directories
install (TARGETS ACTSGenericDetectorExample ACTSGenericDetector
	 EXPORT ACTSTargets
	 LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	 RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install (DIRECTORY include/ACTS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
