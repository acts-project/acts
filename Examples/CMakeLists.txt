# define library target
add_library (ACTSCustomLogger SHARED src/CustomDefaultLogger.cpp)
target_include_directories (ACTSCustomLogger PRIVATE ${CMAKE_SOURCE_DIR}/Core/include)

add_executable (PropagatorExample bin/PropagatorExample.cpp)
target_compile_options (PropagatorExample PRIVATE -g)
target_link_libraries (PropagatorExample PRIVATE ACTSCore)

# set installation directories
install (TARGETS ACTSCustomLogger
         EXPORT ACTSExamplesTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
         COMPONENT Examples)

install (TARGETS PropagatorExample
         RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/Examples)

install (FILES scripts/plot_hits.gp scripts/plot_track.gp
         DESTINATION ${CMAKE_INSTALL_BINDIR}/Examples)

# add to list of supported components
list (APPEND _supported_components Examples)
set (_supported_components "${_supported_components}" PARENT_SCOPE)
