# Configure setup and install python files
set(_python_dir "${CMAKE_BINARY_DIR}/python")
file(MAKE_DIRECTORY ${_python_dir})
file(MAKE_DIRECTORY ${_python_dir}/acts)

set(py_files __init__.py _adapter.py histogram.py)

acts_add_python_module(
    ActsPythonBindings
    SOURCES src/CoreModuleEntry.cpp
    OUTPUT_DIRECTORY ${_python_dir}/acts
    INSTALL_DIRECTORY ${_acts_python_install_dir}
    LIBRARY_DIRECTORY ${CMAKE_INSTALL_LIBDIR}
)

target_sources(
    ActsPythonBindings
    PRIVATE
        src/Definitions.cpp
        src/EventData.cpp
        src/Utilities.cpp
        src/MagneticField.cpp
        src/Material.cpp
        src/Surfaces.cpp
        src/Geometry.cpp
        src/GeometryGen1.cpp
        src/GeometryGen3.cpp
        src/Navigation.cpp
        src/Propagation.cpp
        src/Seeding.cpp
        src/TrackFinding.cpp
        src/Visualization.cpp
)

target_link_libraries(
    ActsPythonBindings
    PUBLIC Acts::Core Acts::PythonUtilities
)

foreach(f ${py_files})
    get_filename_component(_rel ${f} DIRECTORY)
    acts_python_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/python/${f}
        ${_python_dir}/acts/${f}
        ${_acts_python_install_dir}/${_rel}
    )
endforeach()

acts_pybind_generate_stubs(
    acts
    ActsPythonBindings
    PACKAGE
    ENABLED ${ACTS_GENERATE_PYTHON_STUBS}
    PYTHON_ROOT ${_python_dir}
    PYTHON_EXECUTABLE ${Python_EXECUTABLE}
    IGNORE_INVALID_IDENTIFIERS "GeoModelTree::FpvConstLink"
)
