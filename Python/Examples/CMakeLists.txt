set(_python_dir "${CMAKE_BINARY_DIR}/python")
file(MAKE_DIRECTORY ${_python_dir})
file(MAKE_DIRECTORY ${_python_dir}/acts)

set(_python_install_dir "python/acts/examples")
set(py_files
    __init__.py
    simulation.py
    reconstruction.py
    itk.py
    odd.py
)

# Macro to add example bindings conditionally
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

# Add the examples python binding
pybind11_add_module(ActsExamplesPythonBindings src/ExamplesModuleEntry.cpp)

# The PythonBindings for the examples
target_sources(
    ActsExamplesPythonBindings
    PRIVATE
        src/Framework.cpp
        src/ExampleAlgorithms.cpp
        src/Generators.cpp
        src/AmbiguityResolution.cpp
        src/Detector.cpp
        src/Digitization.cpp
        src/MaterialMapping.cpp
        src/Propagation.cpp
        src/TrackFinding.cpp
        src/TrackFitting.cpp
        src/TruthTracking.cpp
        src/Vertexing.cpp
        src/Input.cpp
        src/Output.cpp
        src/Obj.cpp
        src/Utilities.cpp
)

install(TARGETS ActsExamplesPythonBindings DESTINATION ${_python_install_dir})

set_target_properties(
    ActsExamplesPythonBindings
    PROPERTIES INSTALL_RPATH "\$ORIGIN/../../${CMAKE_INSTALL_LIBDIR}"
)
set_target_properties(
    ActsExamplesPythonBindings
    PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_python_dir}/acts/examples
)

target_link_libraries(
    ActsExamplesPythonBindings
    PUBLIC
        Acts::Core
        Acts::ExamplesFramework
        Acts::ExamplesFatras
        Acts::ExamplesAmbiguityResolution
        Acts::ExamplesDetectorsCommon
        Acts::ExamplesDetectorGeneric
        Acts::ExamplesDetectorTelescope
        Acts::ExamplesDigitization
        Acts::ExamplesGenerators
        Acts::ExamplesMaterialMapping
        Acts::ExamplesUtilities
        Acts::ExamplesIoCsv
        Acts::ExamplesIoObj
        Acts::ExamplesPropagation
        Acts::ExamplesTrackFinding
        Acts::ExamplesTrackFitting
        Acts::ExamplesTruthTracking
        Acts::ExamplesVertexing
        Acts::ExamplesPrinters
        Acts::PythonUtilities
)

# Add all plugin related examples bindings conditionally to their build flag
add_examples_binding_if(Alignment "Acts::ExamplesAlignment;Acts::ExamplesDetectorsCommon" ACTS_BUILD_ALIGNMENT FALSE)
add_examples_binding_if(DD4hep Acts::ExamplesDetectorDD4hep ACTS_BUILD_EXAMPLES_DD4HEP TRUE)
add_examples_binding_if(EDM4hep Acts::ExamplesIoEDM4hep ACTS_BUILD_EXAMPLES_EDM4HEP FALSE)
add_examples_binding_if(Geant4
    "Acts::ExamplesGeant4;Acts::ExamplesDetectorGeant4"
    ACTS_BUILD_EXAMPLES_GEANT4
    TRUE
)
add_examples_binding_if(GeoModel
    "Acts::PluginGeoModel;Acts::ExamplesITkModuleSplitting;Acts::ExamplesDetectorGeoModel;Acts::ExamplesMuonSpectrometerMockupDetector"
    ACTS_BUILD_PLUGIN_GEOMODEL TRUE
)
add_examples_binding_if(Gnn Acts::ExamplesTrackFindingGnn ACTS_BUILD_EXAMPLES_GNN FALSE)
add_examples_binding_if(Hashing "Acts::PluginHashing;Acts::ExamplesTrackFinding" ACTS_BUILD_EXAMPLES_HASHING TRUE)
add_examples_binding_if(Json Acts::ExamplesIoJson ACTS_BUILD_PLUGIN_JSON TRUE)
add_examples_binding_if(Onnx "Acts::ExamplesFrameworkML;Acts::ExamplesTrackFindingML" ACTS_BUILD_PLUGIN_ONNX FALSE)
add_examples_binding_if(Pythia8 Acts::ExamplesGeneratorsPythia8 ACTS_BUILD_EXAMPLES_PYTHIA8 FALSE)
add_examples_binding_if(Root Acts::ExamplesIoRoot ACTS_BUILD_PLUGIN_ROOT TRUE)
add_examples_binding_if(Svg Acts::ExamplesIoSvg ACTS_BUILD_PLUGIN_ACTSVG TRUE)
add_examples_binding_if(TGeo Acts::ExamplesDetectorTGeo ACTS_BUILD_PLUGIN_ROOT TRUE)
add_examples_binding_if(Traccc "Acts::PluginCovfie;Acts::ExamplesTraccc" ACTS_BUILD_PLUGIN_TRACCC FALSE)
add_examples_binding_if(TruthJet Acts::ExamplesJets ACTS_BUILD_PLUGIN_FASTJET FALSE)

# Special case for HepMC3 as it a hard dependency
pybind11_add_module(ActsExamplesPythonBindingsHepMC3 src/plugins/HepMC3.cpp)
install(
    TARGETS ActsExamplesPythonBindingsHepMC3 DESTINATION ${_python_install_dir}
)
set_target_properties(
    ActsExamplesPythonBindingsHepMC3
    PROPERTIES INSTALL_RPATH "\$ORIGIN/../../${CMAKE_INSTALL_LIBDIR}"
)
set_target_properties(
    ActsExamplesPythonBindingsHepMC3
    PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_python_dir}/acts/examples
)
target_link_libraries(
    ActsExamplesPythonBindingsHepMC3
    PUBLIC Acts::Core Acts::PythonUtilities Acts::ExamplesIoHepMC3
)
add_dependencies(ActsPythonBindings ActsExamplesPythonBindingsHepMC3)
list(APPEND py_files hepmc3/__init__.py)
list(APPEND py_files hepmc3/__main__.py)

# Symlink hepmc3_normalize.py to build/bin for development
add_custom_target(
    hepmc_normalize_symlink
    ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/bin
    COMMAND
        ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/python/hepmc_normalize.py
        ${CMAKE_BINARY_DIR}/bin/hepmc_normalize.py
    COMMENT "Creating symlink for hepmc_normalize.py"
    VERBATIM
)
add_dependencies(hepmc_normalize_symlink ActsExamplesIoHepMC3)
install(
    PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/python/hepmc_normalize.py
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    RENAME hepmc_normalize.py
)

# Install the additional python files
foreach(f ${py_files})
    set(_target ${_python_dir}/acts/examples/${f})
    get_filename_component(_dir ${_target} DIRECTORY)
    file(MAKE_DIRECTORY ${_dir})

    file(
        CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/python/${f} ${_target}
        SYMBOLIC
    )

    get_filename_component(_rel ${f} DIRECTORY)
    message(
        STATUS
        "Installing example python file: ${f} to ${_python_install_dir}/${_rel}"
    )

    install(
        FILES ${CMAKE_CURRENT_SOURCE_DIR}/python/${f}
        DESTINATION ${_python_install_dir}/${_rel}
    )
endforeach()

# Finally copy the python files
