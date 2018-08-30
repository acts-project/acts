#--- CMake Config Files -----------------------------------------------
# - Use CMake's module to help generating relocatable config files
include (CMakePackageConfigHelpers)

# - Versioning
write_basic_package_version_file ("${CMAKE_CURRENT_BINARY_DIR}/ActsConfigVersion.cmake"
                                  VERSION ${ACTS_VERSION}
                                  COMPATIBILITY SameMajorVersion)

# - Install time config and target files
configure_package_config_file (${CMAKE_CURRENT_LIST_DIR}/ActsConfig.cmake.in "${PROJECT_BINARY_DIR}/ActsConfig.cmake"
                               INSTALL_DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/cmake/Acts"
                               PATH_VARS CMAKE_INSTALL_BINDIR CMAKE_INSTALL_INCLUDEDIR CMAKE_INSTALL_LIBDIR)

# - install and export
install (FILES "${PROJECT_BINARY_DIR}/ActsConfigVersion.cmake" "${PROJECT_BINARY_DIR}/ActsConfig.cmake"
         DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/cmake/Acts")

foreach (_comp ${_supported_components})
  install (EXPORT Acts${_comp}Targets 
           DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/cmake/Acts")
endforeach ()
