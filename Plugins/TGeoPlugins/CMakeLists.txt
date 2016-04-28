# we need root
find_package(ROOT)

# get list of all source files
file(GLOB src_files "src/*.cxx")

# define library target
add_library(ACTSTGeoPlugin SHARED ${src_files})

# setup include directories
target_include_directories(ACTSTGeoPlugin PUBLIC ${ROOT_INCLUDE_DIRS})
target_include_directories(ACTSTGeoPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)

# setup linker dependencies
target_link_libraries(ACTSTGeoPlugin PUBLIC ACTSCore)

# set installation directories
install(TARGETS ACTSTGeoPlugin EXPORT ACTSTargets
	LIBRARY DESTINATION lib)

install(DIRECTORY include/ACTS DESTINATION include)
