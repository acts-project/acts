# Create a setup.sh file at the install prefix to setup the shell environment
configure_file(
  ${CMAKE_CURRENT_LIST_DIR}/setup.sh.in
  ${PROJECT_BINARY_DIR}/setup.sh
  @ONLY)
install(
  FILES ${PROJECT_BINARY_DIR}/setup.sh
  DESTINATION ${CMAKE_INSTALL_PREFIX})
