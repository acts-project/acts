include(FetchContent)
FetchContent_Declare(
    doxygen-awesome-css
    GIT_REPOSITORY https://github.com/jothepro/doxygen-awesome-css.git
    GIT_TAG 1f3620084ff75734ed192101acf40e9dff01d848
)
FetchContent_MakeAvailable(doxygen-awesome-css)

# Save the location the files were cloned into
# This allows us to get the path to doxygen-awesome.css
FetchContent_GetProperties(doxygen-awesome-css SOURCE_DIR AWESOME_CSS_DIR)

set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR})
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

doxygen_add_docs(
    docs
    COMMENT "Generating documentation using doxygen"
    CONFIG_FILE ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
)
