# get list of all source files
file (GLOB_RECURSE src_files "*/*.*pp")

# define library target
add_library (ExamplesLib SHARED ${src_files})
target_include_directories (ExamplesLib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries (ExamplesLib PUBLIC ACore)

# define executables
add_executable (GenericDetector bin/main.cxx)
target_include_directories (GenericDetector PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include/)
target_link_libraries (GenericDetector PRIVATE ACore ExamplesLib)

# set installation directories
install (TARGETS ExamplesLib
         EXPORT ExamplesTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
         COMPONENT ExamplesLib)

install (TARGETS GenericDetector
         RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install (DIRECTORY include/ACTS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
         COMPONENT ExamplesLib)

add_custom_target (ExamplesProject)
add_dependencies (ExamplesProject ExamplesLib GenericDetector)

# add to list of supported components
list (APPEND _supported_components Examples)
set (_supported_components "${_supported_components}" PARENT_SCOPE)
