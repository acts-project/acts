set(doc_sources
  index.rst
  integration/atlas.md
  integration/fcc.md
  integration/integration.rst
  modules/eventdata.md
  modules/geometry.md
  modules/magnetic_field.rst
  modules/material.rst
  modules/modules.rst
  modules/propagation.md
  plugins/dd4hep.rst
  plugins/digitization.rst
  plugins/identification.rst
  plugins/plugins.rst
  plugins/tgeo.rst
  grid_axis.md
  logging.rst
  authors.rst
  ../AUTHORS
  license.rst
  ../LICENSE
  units
)

configure_file(conf.py.in ${CMAKE_CURRENT_BINARY_DIR}/conf.py @ONLY)
configure_file(doxygen.config.in ${CMAKE_CURRENT_BINARY_DIR}/doxygen.config @ONLY)
add_custom_target(
  doc
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
