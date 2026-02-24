unset(IS_INSTALLED_VERSION CACHE)

configure_file(
    ${CMAKE_CURRENT_LIST_DIR}/setup.sh.in
    ${PROJECT_BINARY_DIR}/this_acts.sh
    @ONLY
)
message(STATUS "making setup script with dependencies: ${ROOT_DIR}")
configure_file(
    ${CMAKE_CURRENT_LIST_DIR}/setup_withdeps.sh.in
    ${PROJECT_BINARY_DIR}/this_acts_withdeps.sh
    @ONLY
)

set(IS_INSTALLED_VERSION ON)

configure_file(
    ${CMAKE_CURRENT_LIST_DIR}/setup.sh.in
    ${PROJECT_BINARY_DIR}/install_this_acts.sh
    @ONLY
)
install(
    FILES ${PROJECT_BINARY_DIR}/install_this_acts.sh
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    RENAME this_acts.sh
)

configure_file(
    ${CMAKE_CURRENT_LIST_DIR}/setup_withdeps.sh.in
    ${PROJECT_BINARY_DIR}/install_this_acts_withdeps.sh
    @ONLY
)
install(
    FILES ${PROJECT_BINARY_DIR}/install_this_acts_withdeps.sh
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    RENAME this_acts_withdeps.sh
)

unset(IS_INSTALLED_VERSION CACHE)
