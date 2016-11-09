#--- CMake Config Files -----------------------------------------------
# - Use CMake's module to help generating relocatable config files
include (CMakePackageConfigHelpers)

# - Versioning
write_basic_package_version_file ("${CMAKE_CURRENT_BINARY_DIR}/ACTSConfigVersion.cmake"
                                  VERSION ${ACTS_VERSION}
                                  COMPATIBILITY SameMajorVersion)

# - Install time config and target files
configure_package_config_file (${CMAKE_CURRENT_LIST_DIR}/ACTSConfig.cmake.in "${PROJECT_BINARY_DIR}/ACTSConfig.cmake"
                               INSTALL_DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/cmake/ACTS"
                               PATH_VARS CMAKE_INSTALL_BINDIR CMAKE_INSTALL_INCLUDEDIR CMAKE_INSTALL_LIBDIR)

# - install and export
install (FILES "${PROJECT_BINARY_DIR}/ACTSConfigVersion.cmake" "${PROJECT_BINARY_DIR}/ACTSConfig.cmake"
         DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/cmake/ACTS")

foreach (_comp ${_supported_components})
  install (EXPORT ACTS${_comp}Targets 
           NAMESPACE ACTS:: 
           DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/cmake/ACTS")
endforeach ()

# hack to fix INTERFACE_COMPILE_DEFINITIONS
add_custom_target (fix python ${PROJECT_SOURCE_DIR}/cmake/fix-cmake-target-file.py ACoreTargets.cmake
                   WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATAROOTDIR}/cmake/ACTS VERBATIM)
install (CODE "execute_process(COMMAND \"${CMAKE_COMMAND}\" --build \"${CMAKE_BINARY_DIR}\" --target fix)")
