file(GLOB_RECURSE doxygen_sources
  ${PROJECT_SOURCE_DIR}/Core/include/*.hpp
  ${PROJECT_SOURCE_DIR}/Fatras/include/*.hpp
  ${PROJECT_SOURCE_DIR}/Plugins/*/include/*.hpp
)
set(sphinx_sources
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
  utilities/grid_axis.md
  utilities/logging.rst
  utilities/units.rst
  utilities/utilities.rst
  authors.rst
  codeguide.md
  index.rst
  license.rst
  ../AUTHORS.md
  ../LICENSE
)

set(output_dir ${CMAKE_CURRENT_SOURCE_DIR}/_build)
set(output_html_dir ${output_dir}/html)
set(output_html_index ${output_dir}/html/index.html)

# default target builds the full documentation
add_custom_target(
  doc
  COMMAND
    ${Sphinx_EXECUTABLE}
      -b html
      -d ${output_dir}/doctrees
      -j auto
      -t use_doxygen
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${output_html_dir}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  DEPENDS conf.py Doxyfile ${doxygen_sources} ${sphinx_sources}
  COMMENT "Build documentation")
# helper target for faster iterations without the doxygen API docs
add_custom_target(
  doc-without-doxygen
  COMMAND
    ${Sphinx_EXECUTABLE}
      -b html
      -d ${output_dir}/doctrees
      -j auto
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${output_html_dir}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  DEPENDS conf.py Doxyfile ${doxygen_sources} ${sphinx_sources}
  COMMENT "Build documentation (without Doxygen)")

install(
  DIRECTORY ${output_html_dir}/
  DESTINATION ${CMAKE_INSTALL_DOCDIR}/Acts OPTIONAL)
