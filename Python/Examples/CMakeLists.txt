set(_python_dir "${CMAKE_BINARY_DIR}/python")
file(MAKE_DIRECTORY ${_python_dir})
file(MAKE_DIRECTORY ${_python_dir}/acts)

set(_python_install_dir "python/acts/examples")

macro(
    add_examples_binding_if
    _name
    _link_libraries
    _add_dependencies
    _flag
    _plugin_import
)
    # Create a lower case name version
    string(TOLOWER ${_name} _name_lower)

    if(${_flag})
        # automatically populate the init file with actually built modules
        #list(APPEND _plugins_built "'${_name_lower}'")
        # automatically prefix the target name
        set(_target "ActsExamplesPythonBindings${_name}")
        pybind11_add_module(${_target} src/${_name}.cpp)

        install(TARGETS ${_target} DESTINATION ${_python_install_dir})

        set_target_properties(
            ${_target}
            PROPERTIES INSTALL_RPATH "\$ORIGIN/../../${CMAKE_INSTALL_LIBDIR}"
        )
        set_target_properties(
            ${_target}
            PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_python_dir}/acts/examples
        )

        target_link_libraries(
            ${_target}
            PUBLIC Acts::Core Acts::PythonUtilities ${_link_libraries}
        )
        # Add the dependencies
        add_dependencies(ActsPythonBindings ${_target})
        foreach(_dep IN LISTS ${_add_dependencies})
            add_dependencies(${_dep} ${_target})
        endforeach()
    endif()

    set(EXAMPLE_NAME "${_name}")
    set(EXAMPLE_BUILD_FLAG "${_flag}")
    set(PLUGIN_IMPORT "")
    if(${_plugin_import})
        set(PLUGIN_IMPORT "from acts import ${_name_lower}")
    endif()

    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/python/example.py.in
        ${_python_dir}/acts/examples/${_name_lower}.py
        @ONLY
    )
    install(
        FILES ${_python_dir}/acts/examples/${_name_lower}.py
        DESTINATION ${_python_install_dir}
    )
endmacro()

add_examples_binding_if(DD4hep Acts::ExamplesDetectorDD4hep "" ACTS_BUILD_EXAMPLES_DD4HEP TRUE)
add_examples_binding_if(EDM4hep Acts::ExamplesIoEDM4hep "" ACTS_BUILD_EXAMPLES_EDM4HEP FALSE)
add_examples_binding_if(Geant4
    "Acts::ExamplesGeant4;Acts::ExamplesDetectorGeant4;Acts::ExamplesMuonSpectrometerMockupDetector"
    ""
    ACTS_BUILD_EXAMPLES_GEANT4
    TRUE
)
add_examples_binding_if(Podio Acts::ExamplesIoPodio "" ACTS_BUILD_EXAMPLES_PODIO FALSE)
