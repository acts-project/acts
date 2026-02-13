# Configure setup and install python files
set(_python_dir "${CMAKE_BINARY_DIR}/python")
file(MAKE_DIRECTORY ${_python_dir})
file(MAKE_DIRECTORY ${_python_dir}/acts)

set(py_files __init__.py _adapter.py)

pybind11_add_module(ActsPythonBindings src/CoreModuleEntry.cpp)

target_sources(
    ActsPythonBindings
    PRIVATE
        src/Definitions.cpp
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

install(TARGETS ActsPythonBindings DESTINATION ${ACTS_PYTHON_INSTALL_DIR})

set_target_properties(
    ActsPythonBindings
    PROPERTIES INSTALL_RPATH "\$ORIGIN/../../${CMAKE_INSTALL_LIBDIR}"
)
set_target_properties(
    ActsPythonBindings
    PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_python_dir}/acts
)

target_link_libraries(
    ActsPythonBindings
    PUBLIC Acts::Core Acts::PythonUtilities
)

configure_file(setup.sh.in ${_python_dir}/setup.sh @ONLY)
install(FILES ${_python_dir}/setup.sh DESTINATION python)

foreach(f ${py_files})
    set(_target ${_python_dir}/acts/${f})
    get_filename_component(_dir ${_target} DIRECTORY)
    file(MAKE_DIRECTORY ${_dir})

    file(
        CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/python/${f} ${_target}
        SYMBOLIC
    )

    get_filename_component(_rel ${f} DIRECTORY)

    install(
        FILES ${CMAKE_CURRENT_SOURCE_DIR}/python/${f}
        DESTINATION ${ACTS_PYTHON_INSTALL_DIR}/${_rel}
    )
endforeach()
