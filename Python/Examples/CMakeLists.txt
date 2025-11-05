set(_python_dir "${CMAKE_BINARY_DIR}/python")
file(MAKE_DIRECTORY ${_python_dir})
file(MAKE_DIRECTORY ${_python_dir}/acts)

set(_python_install_dir "python/acts/examples")

macro(add_examples_binding_if _name _link_libraries _flag _plugin_import)
    # Create a lower case name version
    string(TOLOWER ${_name} _name_lower)

    if(${_flag})
        # automatically populate the init file with actually built modules
        #list(APPEND _plugins_built "'${_name_lower}'")
        # automatically prefix the target name
        set(_target "ActsExamplesPythonBindings${_name}")
        pybind11_add_module(${_target} src/plugins/${_name}.cpp)

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
    # Augment the file if needed
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/python/${_name_lower}_aux.py.in")
        file(
            READ
            "${CMAKE_CURRENT_SOURCE_DIR}/python/${_name_lower}_aux.py.in"
            _extra_content
        )
        file(
            APPEND
            "${_python_dir}/acts/examples/${_name_lower}.py"
            "\n${_extra_content}"
        )
    endif()
    # Install the file
    install(
        FILES ${_python_dir}/acts/examples/${_name_lower}.py
        DESTINATION ${_python_install_dir}
    )
endmacro()

add_examples_binding_if(DD4hep Acts::ExamplesDetectorDD4hep ACTS_BUILD_EXAMPLES_DD4HEP TRUE)
add_examples_binding_if(EDM4hep Acts::ExamplesIoEDM4hep ACTS_BUILD_EXAMPLES_EDM4HEP FALSE)
add_examples_binding_if(Geant4
    "Acts::ExamplesGeant4;Acts::ExamplesDetectorGeant4;Acts::ExamplesMuonSpectrometerMockupDetector"
    ACTS_BUILD_EXAMPLES_GEANT4
    TRUE
)
add_examples_binding_if(GeoModel
    "Acts::PluginGeoModel;Acts::ExamplesITkModuleSplitting;Acts::ExamplesDetectorGeoModel;Acts::ExamplesMuonSpectrometerMockupDetector"
    ACTS_BUILD_PLUGIN_GEOMODEL TRUE
)
add_examples_binding_if(Gnn Acts::ExamplesTrackFindingGnn ACTS_BUILD_EXAMPLES_GNN TRUE)
add_examples_binding_if(Hashing "Acts::PluginHashing;Acts::ExamplesTrackFinding" ACTS_BUILD_EXAMPLES_HASHING TRUE)
add_examples_binding_if(HepMC3 Acts::ExamplesIoHepMC3 TRUE FALSE)
add_examples_binding_if(Json Acts::ExamplesIoJson ACTS_BUILD_PLUGIN_JSON TRUE)
add_examples_binding_if(Onnx "Acts::ExamplesFrameworkML;Acts::ExamplesTrackFindingML" ACTS_BUILD_PLUGIN_ONNX FALSE)
add_examples_binding_if(Podio Acts::ExamplesIoPodio ACTS_BUILD_EXAMPLES_PODIO FALSE)
add_examples_binding_if(Pythia8 Acts::ExamplesGeneratorsPythia8 ACTS_BUILD_EXAMPLES_PYTHIA8 FALSE)
add_examples_binding_if(Root Acts::ExamplesIoRoot ACTS_BUILD_PLUGIN_ROOT TRUE)
add_examples_binding_if(Svg Acts::ExamplesIoSvg ACTS_BUILD_PLUGIN_ACTSVG TRUE)
add_examples_binding_if(TGeo Acts::ExamplesDetectorTGeo ACTS_BUILD_PLUGIN_ROOT TRUE)
add_examples_binding_if(Traccc "Acts::PluginCovfie;Acts::ExamplesTraccc" ACTS_BUILD_PLUGIN_TRACCC FALSE)
add_examples_binding_if(TruthJet Acts::ExamplesJets ACTS_BUILD_PLUGIN_FASTJET FALSE)
