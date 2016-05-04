# get list of all source files
file (GLOB src_files "src/*.cxx")

# define library target
add_library(ACTSExamples SHARED ${src_files})

# setup include directories
target_include_directories(ACTSExamples PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)

# set compiler flags
target_compile_options(ACTSExamples PUBLIC -std=c++1y)

# set linker libraries
target_link_libraries(ACTSExamples PUBLIC ACTSCore)

# set installation directories
install (TARGETS ACTSExamples EXPORT ACTSTargets
	 LIBRARY DESTINATION lib
	 RUNTIME DESTINATION bin)

install (DIRECTORY include/ACTS DESTINATION include)
