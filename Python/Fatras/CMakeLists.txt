set(_python_dir "${CMAKE_BINARY_DIR}/python")
file(MAKE_DIRECTORY ${_python_dir})
file(MAKE_DIRECTORY ${_python_dir}/acts)

acts_add_python_module(
    ActsFatrasPythonBindings
    SOURCES src/Fatras.cpp
    OUTPUT_DIRECTORY ${_python_dir}/acts
    INSTALL_DIRECTORY ${_acts_python_install_dir}
    LIBRARY_DIRECTORY ${CMAKE_INSTALL_LIBDIR}
)
target_link_libraries(
    ActsFatrasPythonBindings
    PUBLIC Acts::Core Acts::Fatras Acts::PythonUtilities
)
add_dependencies(ActsPythonBindings ActsFatrasPythonBindings)
acts_pybind_generate_stubs(
    acts.fatras
    ActsFatrasPythonBindings
    ENABLED ${ACTS_GENERATE_PYTHON_STUBS}
    PYTHON_ROOT ${_python_dir}
    PYTHON_EXECUTABLE ${Python_EXECUTABLE}
    DEPENDS ActsPythonBindings
    IGNORE_INVALID_IDENTIFIERS "GeoModelTree::FpvConstLink"
)

acts_python_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/python/fatras.py
    ${_python_dir}/acts/fatras.py
    ${_acts_python_install_dir}
)
