doxygen_add_docs(
    docs
    COMMENT "Generating documentation using doxygen"
    CONFIG_FILE ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
)

include(ExternalProject)

ExternalProject_Add(
    doxygen_awesome_css
    GIT_REPOSITORY https://github.com/jothepro/doxygen-awesome-css.git
    GIT_TAG 9760c30
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/_build/doxygen-awesome-css
)
add_dependencies(docs doxygen_awesome_css)
