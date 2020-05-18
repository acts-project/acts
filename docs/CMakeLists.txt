# to be compatible with readthedocs.org, the documentation build is always
# in-source (as opposed to the configurable out-of-source build directory for
# the code) and is fully driven by the sphinx tool. Doxygen (to extract the
# source code ) is called automatically by sphinx as well. That also means
# that the CMake dependency handling does not apply and is ignored here.

set(output_doctrees ${CMAKE_CURRENT_SOURCE_DIR}/_build/doctrees)
set(output_html ${CMAKE_CURRENT_SOURCE_DIR}/_build/html)

# this is a debug target that is not part of the regular documentation build
# only intended to diagnose issues with the Doxygen configuration
add_custom_target(
  run-doxygen
  COMMAND ${DOXYGEN_EXECUTABLE}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Run Doxygen")

# standard documentation target to build the full documentation
# this can take quite a long time due to the auto-generation of the
add_custom_target(
  docs
  COMMAND
    ${Sphinx_EXECUTABLE}
      -b html
      -d ${output_doctrees}
      -j auto
      -t use_doxygen
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${output_html}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Build documentation")
# helper target for faster iterations without the doxygen API docs
add_custom_target(
  docs-without-doxygen
  COMMAND
    ${Sphinx_EXECUTABLE}
      -b html
      -d ${output_doctrees}
      -j auto
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${output_html}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Build documentation (without Doxygen)")

install(
  DIRECTORY ${output_html}/
  DESTINATION ${CMAKE_INSTALL_DOCDIR}/Acts OPTIONAL)
