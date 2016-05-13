# get list of all source files
file (GLOB src_files "src/*.cpp")

# define library target
add_library (ACTSGenericDetectorPlugin SHARED ${src_files})

# setup include directories
target_include_directories (ACTSGenericDetectorPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)

# setup linker dependencies
target_link_libraries (ACTSGenericDetectorPlugin PUBLIC ACTSCore)

# set installation directories
install (TARGETS ACTSGenericDetectorPlugin EXPORT ACTSTargets
	 LIBRARY DESTINATION lib)

install (DIRECTORY include/ACTS DESTINATION include)
