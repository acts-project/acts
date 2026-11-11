set(_python_dir "${CMAKE_BINARY_DIR}/python")
file(MAKE_DIRECTORY ${_python_dir})
file(MAKE_DIRECTORY ${_python_dir}/acts)

set(_examples_python_install_dir "${_acts_python_install_dir}/examples")
set(py_files
    __init__.py
    simulation.py
    reconstruction.py
    itk.py
    odd.py
    dataset.py
    uproot.py
    visualization.py
    scipy.py
)

# Macro to add example bindings conditionally
macro(add_examples_binding_if)
    set(_options)
    set(_one_value_args NAME FLAG PLUGIN_IMPORT)
    set(_multi_value_args LINK_LIBRARIES)
    cmake_parse_arguments(
        AEBI
        "${_options}"
        "${_one_value_args}"
        "${_multi_value_args}"
        ${ARGN}
    )

    if(AEBI_UNPARSED_ARGUMENTS)
        message(
            FATAL_ERROR
            "add_examples_binding_if got unexpected arguments: ${AEBI_UNPARSED_ARGUMENTS}"
        )
    endif()

    if(
        NOT DEFINED AEBI_NAME
        OR NOT DEFINED AEBI_FLAG
        OR NOT DEFINED AEBI_PLUGIN_IMPORT
    )
        message(
            FATAL_ERROR
            "add_examples_binding_if requires NAME, FLAG, and PLUGIN_IMPORT arguments."
        )
    endif()

    # Create a lower case name version
    string(TOLOWER ${AEBI_NAME} _name_lower)

    if(${AEBI_FLAG})
        # automatically populate the init file with actually built modules
        #list(APPEND _plugins_built "'${_name_lower}'")
        # automatically prefix the target name
        set(_target "ActsExamplesPythonBindings${AEBI_NAME}")

        acts_add_python_module(
            ${_target}
            SOURCES src/plugins/${AEBI_NAME}.cpp
            OUTPUT_DIRECTORY ${_python_dir}/acts/examples
            INSTALL_DIRECTORY ${_examples_python_install_dir}
            LIBRARY_DIRECTORY ${CMAKE_INSTALL_LIBDIR}
        )

        target_link_libraries(
            ${_target}
            PUBLIC Acts::Core Acts::PythonUtilities ${AEBI_LINK_LIBRARIES}
        )
        # Add the dependencies
        add_dependencies(ActsPythonBindings ${_target})
    endif()

    set(EXAMPLE_NAME "${AEBI_NAME}")
    set(EXAMPLE_BUILD_FLAG "${AEBI_FLAG}")
    set(PLUGIN_IMPORT "")
    if(${AEBI_PLUGIN_IMPORT})
        set(PLUGIN_IMPORT "from acts import ${_name_lower}")
    endif()

    acts_python_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/python/example.py.in
        ${_python_dir}/acts/examples/${_name_lower}.py
        ${_examples_python_install_dir}
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
endmacro()

# Add the examples python binding
acts_add_python_module(
    ActsExamplesPythonBindings
    SOURCES src/ExamplesModuleEntry.cpp
    OUTPUT_DIRECTORY ${_python_dir}/acts/examples
    INSTALL_DIRECTORY ${_examples_python_install_dir}
    LIBRARY_DIRECTORY ${CMAKE_INSTALL_LIBDIR}
)

add_dependencies(ActsPythonBindings ActsExamplesPythonBindings)

# The PythonBindings for the examples
target_sources(
    ActsExamplesPythonBindings
    PRIVATE
        src/Framework.cpp
        src/EventData.cpp
        src/EventDataMeasurements.cpp
        src/EventDataTrackContainers.cpp
        src/EventDataTracks.cpp
        src/EventDataTrackStates.cpp
        src/ExampleAlgorithms.cpp
        src/Generators.cpp
        src/AmbiguityResolution.cpp
        src/Detector.cpp
        src/Digitization.cpp
        src/MaterialMapping.cpp
        src/Propagation.cpp
        src/PythonSpecific.cpp
        src/TrackFinding.cpp
        src/TrackFitting.cpp
        src/TruthTracking.cpp
        src/HelloWorld.cpp
        src/Vertexing.cpp
        src/Input.cpp
        src/Output.cpp
        src/Obj.cpp
        src/Utilities.cpp
)

acts_pybind_generate_stubs(
    acts.examples
    ActsExamplesPythonBindings
    PACKAGE
    ENABLED ${ACTS_GENERATE_PYTHON_STUBS}
    PYTHON_ROOT ${_python_dir}
    PYTHON_EXECUTABLE ${Python_EXECUTABLE}
    DEPENDS ActsPythonBindings
    IGNORE_INVALID_IDENTIFIERS "GeoModelTree::FpvConstLink"
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
        Acts::ExamplesHelloWorld
        Acts::ExamplesVertexing
        Acts::ExamplesPrinters
        Acts::PythonUtilities
)

# Add all plugin related examples bindings conditionally to their build flag
add_examples_binding_if(
    NAME Alignment
    LINK_LIBRARIES Acts::ExamplesAlignment Acts::ExamplesDetectorsCommon
    FLAG ACTS_BUILD_ALIGNMENT
    PLUGIN_IMPORT FALSE
)
add_examples_binding_if(
    NAME AlignmentMillePede
    LINK_LIBRARIES Acts::ExamplesAlignmentMillePede
    FLAG ACTS_BUILD_PLUGIN_MILLE
    PLUGIN_IMPORT FALSE
)
add_examples_binding_if(
    NAME Arrow
    LINK_LIBRARIES Acts::ExamplesIoParquet Acts::ExamplesIoArrow
    FLAG ACTS_BUILD_EXAMPLES_PARQUET
    PLUGIN_IMPORT TRUE
)
add_examples_binding_if(
    NAME DD4hep
    LINK_LIBRARIES Acts::ExamplesDetectorDD4hep
    FLAG ACTS_BUILD_EXAMPLES_DD4HEP
    PLUGIN_IMPORT TRUE
)
add_examples_binding_if(
    NAME EDM4hep
    LINK_LIBRARIES Acts::ExamplesIoEDM4hep
    FLAG ACTS_BUILD_EXAMPLES_EDM4HEP
    PLUGIN_IMPORT FALSE
)
add_examples_binding_if(
    NAME Geant4
    LINK_LIBRARIES Acts::ExamplesGeant4 Acts::ExamplesDetectorGeant4
    FLAG ACTS_BUILD_EXAMPLES_GEANT4
    PLUGIN_IMPORT TRUE
)
add_examples_binding_if(
    NAME GeoModel
    LINK_LIBRARIES
        Acts::PluginGeoModel
        Acts::ExamplesITkModuleSplitting
        Acts::ExamplesDetectorGeoModel
        Acts::ExamplesMuonSpectrometerMockupDetector
    FLAG ACTS_BUILD_PLUGIN_GEOMODEL
    PLUGIN_IMPORT TRUE
)
add_examples_binding_if(
    NAME Gnn
    LINK_LIBRARIES Acts::ExamplesTrackFindingGnn
    FLAG ACTS_BUILD_EXAMPLES_GNN
    PLUGIN_IMPORT TRUE
)
add_examples_binding_if(
    NAME Hashing
    LINK_LIBRARIES Acts::ExamplesTrackFinding
    FLAG ACTS_BUILD_EXAMPLES_HASHING
    PLUGIN_IMPORT FALSE
)
add_examples_binding_if(
    NAME Json
    LINK_LIBRARIES Acts::ExamplesIoJson
    FLAG ACTS_BUILD_PLUGIN_JSON
    PLUGIN_IMPORT TRUE
)
add_examples_binding_if(
    NAME Onnx
    LINK_LIBRARIES Acts::ExamplesTrackFindingML
    FLAG ACTS_BUILD_PLUGIN_ONNX
    PLUGIN_IMPORT FALSE
)
add_examples_binding_if(
    NAME Pythia8
    LINK_LIBRARIES Acts::ExamplesGeneratorsPythia8
    FLAG ACTS_BUILD_EXAMPLES_PYTHIA8
    PLUGIN_IMPORT FALSE
)
add_examples_binding_if(
    NAME Root
    LINK_LIBRARIES Acts::ExamplesIoRoot
    FLAG ACTS_BUILD_EXAMPLES_ROOT
    PLUGIN_IMPORT TRUE
)
add_examples_binding_if(
    NAME Svg
    LINK_LIBRARIES Acts::ExamplesIoSvg
    FLAG ACTS_BUILD_PLUGIN_ACTSVG
    PLUGIN_IMPORT TRUE
)
add_examples_binding_if(
    NAME TGeo
    LINK_LIBRARIES Acts::ExamplesDetectorTGeo
    FLAG ACTS_BUILD_EXAMPLES_ROOT
    PLUGIN_IMPORT TRUE
)
add_examples_binding_if(
    NAME Detray
    LINK_LIBRARIES Acts::PluginCovfie Acts::ExamplesDetray
    FLAG ACTS_BUILD_PLUGIN_TRACCC
    PLUGIN_IMPORT FALSE
)
add_examples_binding_if(
    NAME TruthJet
    LINK_LIBRARIES Acts::ExamplesJets
    FLAG ACTS_BUILD_EXAMPLES_FASTJET
    PLUGIN_IMPORT FALSE
)

# Special case for HepMC3 as it a hard dependency
acts_add_python_module(
    ActsExamplesPythonBindingsHepMC3
    SOURCES src/plugins/HepMC3.cpp
    OUTPUT_DIRECTORY ${_python_dir}/acts/examples
    INSTALL_DIRECTORY ${_examples_python_install_dir}
    LIBRARY_DIRECTORY ${CMAKE_INSTALL_LIBDIR}
)
target_link_libraries(
    ActsExamplesPythonBindingsHepMC3
    PUBLIC Acts::Core Acts::PythonUtilities Acts::ExamplesIoHepMC3
)
add_dependencies(ActsPythonBindings ActsExamplesPythonBindingsHepMC3)
list(APPEND py_files hepmc3/__init__.py)
list(APPEND py_files hepmc3/__main__.py)

# Handle python scripts
file(GLOB _python_scripts ${CMAKE_CURRENT_SOURCE_DIR}/bin/*.py)

# Symlink scripts to build dir. Create a custom command for each script
# and a custom target that depends on all of them.
set(_python_scripts_symlinks)
foreach(_src ${_python_scripts})
    get_filename_component(_name ${_src} NAME)
    set(_dst ${CMAKE_BINARY_DIR}/bin/${_name})
    add_custom_command(
        OUTPUT ${_dst}
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/bin
        COMMAND ${CMAKE_COMMAND} -E create_symlink ${_src} ${_dst}
        DEPENDS ${_src}
        COMMENT "Creating symlink for ${_script}"
        VERBATIM
    )
    list(APPEND _python_scripts_symlinks ${_dst})
endforeach()

add_custom_target(python_bin_symlink ALL DEPENDS ${_python_scripts_symlinks})
add_dependencies(python_bin_symlink ActsExamplesPythonBindings)

# Install scripts from bin directory
install(PROGRAMS ${_python_scripts} DESTINATION ${CMAKE_INSTALL_BINDIR})

# Install the additional python files
foreach(f ${py_files})
    get_filename_component(_rel ${f} DIRECTORY)
    acts_python_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/python/${f}
        ${_python_dir}/acts/examples/${f}
        ${_examples_python_install_dir}/${_rel}
    )
endforeach()

# Finally copy the python files
