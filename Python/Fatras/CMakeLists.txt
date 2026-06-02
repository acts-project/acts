set(_python_dir "${CMAKE_BINARY_DIR}/python")
file(MAKE_DIRECTORY ${_python_dir})
file(MAKE_DIRECTORY ${_python_dir}/acts)

pybind11_add_module(ActsFatrasPythonBindings src/Fatras.cpp)
install(
    TARGETS ActsFatrasPythonBindings
    DESTINATION ${_acts_python_install_dir}
)

set_target_properties(
    ActsFatrasPythonBindings
    PROPERTIES INSTALL_RPATH "${_acts_python_modules_rpath}"
)
set_target_properties(
    ActsFatrasPythonBindings
    PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_python_dir}/acts
)

target_link_libraries(
    ActsFatrasPythonBindings
    PUBLIC Acts::Core Acts::Fatras Acts::PythonUtilities
)
add_dependencies(ActsPythonBindings ActsFatrasPythonBindings)
acts_pybind_generate_stubs(acts.fatras ActsFatrasPythonBindings ActsPythonBindings)

set(_target ${_python_dir}/acts/fatras.py)
file(
    CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/python/fatras.py ${_target}
    SYMBOLIC
)
install(
    FILES ${CMAKE_CURRENT_SOURCE_DIR}/python/fatras.py
    DESTINATION ${_acts_python_install_dir}
)

set(py_files fatras.py)

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
        DESTINATION ${_acts_python_install_dir}/${_rel}
    )
endforeach()
