# to be compatible with readthedocs.org, the documentation build is always
# in-source; as opposed to the configurable out-of-source build directory for
# the code. when running running on readthedocs.org, the build is fully driven
# by Sphinx, including running Doxygen.
#
# this CMake-based build is only intended for local development. Doxygen is
# run separate from Sphinx to avoid running it when no code has change. since
# the full documentation with automatic generation of API docs take O(10mins)
# to build, the default `docs` target does not include the API docs generation
# but includes the Sphinx-Doxygen integration (breathe).
#
# the full documentation as on readthedocs.org can be build via the additional
# `docs-with-api` target.
#
# WARNING: after you ran `docs-with-api` the regular `docs` target will
#          accidentaly pick up the auto-generated api files. remove the
#          `_build` and the `api` directory manually to avoid this.

# this should match the INPUT declaration in the Doxyfile
# CONFIGURE_DEPENDS ensures that adding/removing files leads to reprocessing
file(
  GLOB_RECURSE doxygen_sources
  ../Core/include/*.hpp
  ../Fatras/include/*.hpp
  ../Plugins/include/*.hpp)

set(doxygen_index ${CMAKE_CURRENT_SOURCE_DIR}/_build/doxygen-xml/index.xml)
set(sphinx_doctrees ${CMAKE_CURRENT_SOURCE_DIR}/_build/doctrees)
set(sphinx_html ${CMAKE_CURRENT_SOURCE_DIR}/_build/html)

add_custom_command(
  OUTPUT ${doxygen_index}
  COMMAND ${DOXYGEN_EXECUTABLE}
  DEPENDS ${doxygen_sources} Doxyfile
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Parse source code documentation with Doxygen")

# this is a debug target to diagnose issues with the Doxygen configuration
add_custom_target(run-doxygen DEPENDS ${doxygen_index})

# standard target to build the documentation without automatic API generation
add_custom_target(
  docs
  COMMAND
    ${Sphinx_EXECUTABLE}
      -b html
      -d ${sphinx_doctrees}
      -j auto
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${sphinx_html}
  DEPENDS ${doxygen_index}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Build documentation WITHOUT API documentation")
# extended target to build the full documentation with automatic API generation
add_custom_target(
  docs-with-api
  COMMAND
    ${Sphinx_EXECUTABLE}
      -b html
      -d ${sphinx_doctrees}
      -j auto
      -t use_exhale
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${sphinx_html}
  DEPENDS ${doxygen_index}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Build full documentation")

install(
  DIRECTORY ${sphinx_html}/
  DESTINATION ${CMAKE_INSTALL_DOCDIR}/Acts OPTIONAL)
