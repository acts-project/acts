# get list of all source files
file (GLOB src_files "src/*.cxx")

# define library target
add_library(ACTSGenericDetectorExample SHARED ${src_files})

# setup include directories
target_include_directories(ACTSGenericDetectorExample PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)

# setup linker dependencies
target_link_libraries(ACTSGenericDetectorExample PUBLIC ACTSCore)
target_link_libraries(ACTSGenericDetectorExample PUBLIC ACTSGenericDetectorPlugins)

# set installation directories
install(TARGETS ACTSGenericDetectorExample EXPORT ACTSTargets
LIBRARY DESTINATION lib)

install(DIRECTORY include/ACTS DESTINATION include)
