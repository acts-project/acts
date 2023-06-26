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
#          accidentally pick up the auto-generated api files. remove the
#          `_build` and the `api` directory manually to avoid this.

# this should match the INPUT declaration in the Doxyfile
# CONFIGURE_DEPENDS ensures that adding/removing files leads to reprocessing
file(
  GLOB_RECURSE doxygen_sources
  ../Core/include/*.hpp
  ../Fatras/include/*.hpp
  ../Plugins/include/*.hpp)

set(sphinx_build ${CMAKE_CURRENT_SOURCE_DIR}/_build)
set(doxygen_index ${CMAKE_CURRENT_SOURCE_DIR}/_build/doxygen-xml/index.xml)
set(sphinx_doctrees ${CMAKE_CURRENT_SOURCE_DIR}/_build/doctrees)
set(sphinx_html ${CMAKE_CURRENT_SOURCE_DIR}/_build/html)
set(sphinx_api ${CMAKE_CURRENT_SOURCE_DIR}/api)

add_custom_command(
  OUTPUT ${doxygen_index}
  COMMAND ${DOXYGEN_EXECUTABLE}
  DEPENDS ${doxygen_sources} Doxyfile
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Parse source code documentation with Doxygen")

# this is a debug target to diagnose issues with the Doxygen configuration
add_custom_target(run-doxygen DEPENDS ${doxygen_index})

set(api_dir ${CMAKE_CURRENT_SOURCE_DIR}/api)
set(api_index ${api_dir}/api.rst)

set(doxygen_command)
# list(APPEND doxygen_command "${Sphinx_EXECUTABLE}")
list(APPEND doxygen_command "-b;html")
list(APPEND doxygen_command "-d;${sphinx_doctrees}")
list(APPEND doxygen_command "-j;auto")
list(APPEND doxygen_command "-W;--keep-going")

if(DEFINED ENV{CI})
  list(APPEND doxygen_command "-b;linkcheck")
endif()

list(APPEND doxygen_command "${CMAKE_CURRENT_SOURCE_DIR}")
list(APPEND doxygen_command "${sphinx_html}")

# standard target to build the documentation without automatic API generation
add_custom_target(
  docs
  COMMAND ${Sphinx_EXECUTABLE} ${doxygen_command}
  DEPENDS ${doxygen_index}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Build documentation WITHOUT API documentation")

list(APPEND doxygen_command "-t;run_apidoc")

# extended target to build the full documentation with automatic API generation
add_custom_target(
  docs-with-api
  COMMAND ${Sphinx_EXECUTABLE} ${doxygen_command}
  DEPENDS ${doxygen_index}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Build full documentation")

add_custom_target(
  docs-clean
  COMMAND find ${sphinx_api} -mindepth 1 -type f -and -not -name "api_stub.rst" -and -not -name "api_index.rst" -delete
  COMMAND find ${sphinx_api} -mindepth 1 -type d -delete
  COMMAND rm -rf ${sphinx_build}
  COMMENT "Cleaning documentation artifacts"
)

install(
  DIRECTORY ${sphinx_html}/
  DESTINATION ${CMAKE_INSTALL_DOCDIR}/Acts OPTIONAL)
