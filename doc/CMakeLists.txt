set(doc_sources
  index.rst
  annotated_api.rst
  grid_axis.md)

configure_file(conf.py.in ${CMAKE_CURRENT_BINARY_DIR}/conf.py @ONLY)
configure_file(doxygen.config.in ${CMAKE_CURRENT_BINARY_DIR}/doxygen.config @ONLY)
# manual mkdir since doxygen can not created nested directories. wtf?
add_custom_target(
  doc
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/html/doxygen
  COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doxygen.config
  COMMAND ${Sphinx_EXECUTABLE} -b html -d ${CMAKE_CURRENT_BINARY_DIR}/doctrees -j auto -c ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/html
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generate html documentation with Doxygen and Sphinx"
  VERBATIM
  SOURCES
    doxygen_layout.xml
    ${CMAKE_CURRENT_BINARY_DIR}/doxygen.config
    ${CMAKE_CURRENT_BINARY_DIR}/conf.py
    ${doc_sources})

install(
  DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html
  DESTINATION ${CMAKE_INSTALL_DOCDIR} OPTIONAL)
