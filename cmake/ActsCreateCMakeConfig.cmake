# Create CMake relocatable config files

include(CMakePackageConfigHelpers)

# version is taken automatically from PROJECT_VERSION; no need to specify
write_basic_package_version_file(
  ${PROJECT_BINARY_DIR}/ActsConfigVersion.cmake
  COMPATIBILITY SameMajorVersion)
configure_package_config_file(
  ${CMAKE_CURRENT_LIST_DIR}/ActsConfig.cmake.in
  ${PROJECT_BINARY_DIR}/ActsConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/Acts
  PATH_VARS CMAKE_INSTALL_BINDIR CMAKE_INSTALL_INCLUDEDIR CMAKE_INSTALL_LIBDIR)

# install core cmake configs
install(
  FILES
    ${PROJECT_BINARY_DIR}/ActsConfigVersion.cmake
    ${PROJECT_BINARY_DIR}/ActsConfig.cmake
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/Acts)
# install target configs for all available components
foreach(_component ${_supported_components})
  install(
    EXPORT Acts${_component}Targets
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/Acts)
endforeach()
