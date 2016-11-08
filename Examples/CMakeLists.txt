# get list of all source files
file (GLOB_RECURSE src_files "*/*.*pp")

# define library target
add_library (ACTSExamples SHARED ${src_files})
target_include_directories (ACTSExamples PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries (ACTSExamples PUBLIC ACTSCore)

# define executables
add_executable (ACTSGenericDetector bin/main.cxx)
target_include_directories (ACTSGenericDetector PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include/)
target_link_libraries (ACTSGenericDetector PRIVATE ACTSCore ACTSExamples)

# set installation directories
install (TARGETS ACTSExamples
	 EXPORT ACTSExamplesTargets
	 LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	 COMPONENT Examples)

install (TARGETS ACTSGenericDetector
	 RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install (DIRECTORY include/ACTS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
	 COMPONENT Examples)

# add to list of supported components
list (APPEND _supported_components Examples)
set (_supported_components "${_supported_components}" PARENT_SCOPE)

# add dummy project capturing all targets for the CDash Examples subproject
add_custom_target (ExamplesProject)
add_dependencies (ExamplesProject ACTSExamples ACTSGenericDetector)
