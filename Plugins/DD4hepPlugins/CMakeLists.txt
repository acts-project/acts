# we need DD4Hep
find_package(DD4hep REQUIRED DDCore DDSegmentation)

# get list of all source files
file (GLOB src_files "src/*.cxx")

# define library target
add_library(ACTSDD4hepPlugin SHARED ${src_files})

# setup include directories
target_include_directories(ACTSDD4hepPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_include_directories(ACTSDD4hepPlugin PUBLIC ${DD4hep_INCLUDE_DIRS})

# setup linker dependencies
target_link_libraries(ACTSDD4hepPlugin PUBLIC ACTSCore)
target_link_libraries(ACTSDD4hepPlugin PUBLIC ACTSTGeoPlugin)

# set installation directories
install(TARGETS ACTSDD4hepPlugin EXPORT ACTSTargets
	LIBRARY DESTINATION lib)

install(DIRECTORY include/ACTS DESTINATION include)

