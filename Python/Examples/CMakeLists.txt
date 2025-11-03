set(_python_dir "${CMAKE_BINARY_DIR}/python")
file(MAKE_DIRECTORY ${_python_dir})
file(MAKE_DIRECTORY ${_python_dir}/acts)

set(_python_install_dir "python/acts")

# DD4hep
if(ACTS_BUILD_PLUGIN_DD4HEP AND ACTS_BUILD_EXAMPLES_DD4HEP)
    pybind11_add_module(ActsExamplesPythonBindingsDD4hep src/DD4hepComponent.cpp)
    target_link_libraries(
        ActsExamplesPythonBindingsDD4hep
        PUBLIC Acts::PythonUtilities Acts::ExamplesDetectorDD4hep
    )
    add_dependencies(ActsPythonBindings ActsExamplesPythonBindingsDD4hep)

    install(
        TARGETS
            ActsExamplesPythonBindingsDD4hep
            DESTINATION
            ${_python_install_dir}
    )
    set_target_properties(
        ActsExamplesPythonBindingsDD4hep
        PROPERTIES INSTALL_RPATH "\$ORIGIN/../../${CMAKE_INSTALL_LIBDIR}"
    )
    set_target_properties(
        ActsExamplesPythonBindingsDD4hep
        PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_python_dir}/acts
    )
    list(APPEND py_files examples/dd4hep.py)
endif()

# EDM4hep
if(ACTS_BUILD_EXAMPLES_EDM4HEP)
    pybind11_add_module(ActsExamplesPythonBindingsEDM4hep src/EDM4hepComponent.cpp)
    target_link_libraries(
        ActsExamplesPythonBindingsEDM4hep
        PUBLIC Acts::PythonUtilities Acts::ExamplesIoEDM4hep
    )
    add_dependencies(ActsPythonBindings ActsExamplesPythonBindingsEDM4hep)
    install(
        TARGETS ActsExamplesPythonBindingsEDM4hep DESTINATION ${_python_install_dir}
    )
    set_target_properties(
        ActsExamplesPythonBindingsEDM4hep
        PROPERTIES INSTALL_RPATH "\$ORIGIN/../../${CMAKE_INSTALL_LIBDIR}"
    )
    set_target_properties(
        ActsExamplesPythonBindingsEDM4hep
        PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_python_dir}/acts
    )
    list(APPEND py_files examples/edm4hep.py)
endif()

# Geant4
if(ACTS_BUILD_EXAMPLES_GEANT4)
    pybind11_add_module(ActsExamplesPythonBindingsGeant4 src/Geant4Component.cpp)
    target_link_libraries(
        ActsExamplesPythonBindingsGeant4
        PUBLIC
            Acts::ExamplesGeant4
            Acts::ExamplesDetectorGeant4
            Acts::PythonUtilities
            Acts::ExamplesMuonSpectrometerMockupDetector
    )
    add_dependencies(ActsPythonBindings ActsExamplesPythonBindingsGeant4)

    install(
        TARGETS
            ActsExamplesPythonBindingsGeant4
            DESTINATION
            ${_python_install_dir}
    )
    set_target_properties(
        ActsExamplesPythonBindingsGeant4
        PROPERTIES INSTALL_RPATH "\$ORIGIN/../../${CMAKE_INSTALL_LIBDIR}"
    )
    set_target_properties(
        ActsExamplesPythonBindingsGeant4
        PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_python_dir}/acts
    )
    list(APPEND py_files examples/geant4/__init__.py)
endif()


# Podio
if(ACTS_BUILD_EXAMPLES_PODIO)
    pybind11_add_module(ActsPythonBindingsPodio src/PodioComponent.cpp)
    target_link_libraries(
        ActsPythonBindingsPodio
        PUBLIC Acts::PythonUtilities Acts::ExamplesIoPodio
    )
    add_dependencies(ActsPythonBindings ActsPythonBindingsPodio)
    install(TARGETS ActsPythonBindingsPodio DESTINATION ${_python_install_dir})
    set_target_properties(
        ActsPythonBindingsPodio
        PROPERTIES INSTALL_RPATH "\$ORIGIN/../../${CMAKE_INSTALL_LIBDIR}"
    )
    set_target_properties(
        ActsPythonBindingsPodio
        PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_python_dir}/acts
    )
    list(APPEND py_files examples/podio.py)
endif()
