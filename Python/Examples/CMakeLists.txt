if(ACTS_BUILD_EXAMPLES AND FALSE)
    # Add all components that do not need extra dependencies
    target_sources(
        ActsPythonBindings
        PRIVATE
        src/ExamplesModuleEntry.cpp
        src/Framework.cpp
        src/Detector.cpp
        src/Generators.cpp
        src/MagneticFieldMaps.cpp
        src/AmbiguityResolution.cpp
        src/DigitizationAlgorithms.cpp
        src/MaterialMappingAlgorithms.cpp
        src/TrackFindingAlgorithms.cpp
        src/TrackFittingAlgorithms.cpp
        src/TruthTrackingAlgorithms.cpp
        src/VertexingAlgorithms.cpp
        src/PrinterAlgorithms.cpp
        src/PropagationAlgorithms.cpp
        src/UtilityAlgorithms.cpp
        src/Input.cpp
        src/Output.cpp
    )

    target_link_libraries(
        ActsPythonBindings
        PUBLIC
        ActsExamplesFramework
        ActsExamplesGenerators
        ActsExamplesDetectorsCommon
        ActsExamplesDetectorGeneric
        ActsExamplesDetectorTelescope
        ActsExamplesAmbiguityResolution
        ActsExamplesDigitization
        ActsExamplesMagneticField
        ActsExamplesMaterialMapping
        ActsExamplesTrackFinding
        ActsExamplesTrackFitting
        ActsExamplesTruthTracking
        ActsExamplesVertexing
        ActsExamplesPrinters
        ActsExamplesPropagation
        ActsExamplesUtilities
        ActsExamplesIoCsv
        ActsExamplesIoObj
    )

    # Plugin dependent components, attached to the ActsPythonBindings target

    # Fatras
    #if(ACTS_BUILD_PLUGIN_FATRAS)
    #    target_link_libraries(ActsPythonBindings PUBLIC ActsExamplesFatras)
    #    target_sources(ActsPythonBindings PRIVATE src/FatrasAlgorithms.cpp)
    #else()
    #    target_sources(ActsPythonBindings PRIVATE stubs/FatrasAlgorithmsStub.cpp)
    #endif()

    # GeoModel
    #if(ACTS_BUILD_PLUGIN_GEOMODEL)
    #    target_link_libraries(
    #        ActsPythonBindings PUBLIC
    #        ActsExamplesDetectorGeoModel
    #        ActsExamplesMuonSpectrometerMockupDetector
    #        ActsExamplesITkModuleSplitting
    #    )
    #    target_sources(
    #        ActsPythonBindings PRIVATE src/GeoModelDetector.cpp
    #    )
    #else()
    #    target_sources(
    #        ActsPythonBindings PRIVATE stubs/GeoModelDetectorStub.cpp
    #    )
    #endif()

    # ROOT
    #if(ACTS_BUILD_PLUGIN_ROOT)
    #    target_link_libraries(ActsPythonBindings PUBLIC ActsPluginRoot ActsExamplesDetectorTGeo ActsExamplesIoRoot)
    #    target_sources(ActsPythonBindings PRIVATE src/TGeoDetector.cpp src/RootInput.cpp src/RootOutput.cpp)
    #else()
    #    target_sources(ActsPythonBindings PRIVATE stubs/TGeoDetectorStub.cpp stubs/RootStub.cpp)
    #endif()

    # JSON
    #if(ACTS_BUILD_PLUGIN_JSON)
    #    target_link_libraries(ActsPythonBindings PUBLIC ActsExamplesIoJson)
    #    target_sources(ActsPythonBindings PRIVATE src/JsonInputOutput.cpp)
    #else()
    #    target_sources(ActsPythonBindings PRIVATE stubs/JsonInputOutputStub.cpp)
    #endif()

    # SVG
    #if(ACTS_BUILD_PLUGIN_ACTSVG)
    #    target_link_libraries(ActsPythonBindings PUBLIC ActsExamplesIoSvg)
    #    target_sources(ActsPythonBindings PRIVATE src/SvgOutput.cpp)
    #else()
    #    target_sources(ActsPythonBindings PRIVATE stubs/SvgOutputStub.cpp)
    #endif()

    # traccc
    #if(ACTS_BUILD_PLUGIN_TRACCC)
    #    target_link_libraries(ActsPythonBindings PUBLIC ActsExamplesTraccc)
    #    target_sources(ActsPythonBindings PRIVATE src/TracccAlgorithms.cpp)
    #else()
    #    target_sources(ActsPythonBindings PRIVATE stubs/TracccAlgorithmsStub.cpp)
    #endif()

    # Pythia8
    #if(ACTS_BUILD_PLUGIN_PYTHIA8)
    #    target_link_libraries(ActsPythonBindings PUBLIC ActsExamplesPythia8)
    #    target_sources(ActsPythonBindings PRIVATE src/Pythia8Generator.cpp)
    #else()
    #    target_sources(ActsPythonBindings PRIVATE stubs/Pythia8GeneratorStub.cpp)
    #endif()

    # ExaTrk.X
    #if(ACTS_BUILD_PLUGIN_EXATRKX)
    #    target_link_libraries(ActsPythonBindings PUBLIC ActsExamplesExaTrkX)
    #    target_sources(ActsPythonBindings PRIVATE src/ExaTrkXAlgorithms.cpp)
    #else()
    #    target_sources(ActsPythonBindings PRIVATE stubs/ExaTrkXAlgorithmsStub.cpp)
    #endif()

    # Onnx
    #if(ACTS_BUILD_PLUGIN_ONNX)
    #    target_link_libraries(ActsPythonBindings PUBLIC ActsExamplesOnnx)
    #    target_sources(ActsPythonBindings PRIVATE src/OnnxAlgorithms.cpp)
    #else()
    #    target_sources(ActsPythonBindings PRIVATE stubs/OnnxAlgorithmsStub.cpp)
    #endif()

    # TruthJet
    #if(ACTS_BUILD_PLUGIN_FASTJET)
    #    target_link_libraries(ActsPythonBindings PUBLIC ActsExamplesJets)
    #    target_sources(ActsPythonBindings PRIVATE src/TruthJetAlgorithms.cpp)
    #else()
    #    target_sources(ActsPythonBindings PRIVATE stubs/TruthJetAlgorithmsStub.cpp)
    #endif()

    #if(ACTS_BUILD_EXAMPLES_HASHING)
    #    target_link_libraries(ActsPythonBindings PUBLIC ActsPluginHashing)
    #    target_sources(ActsPythonBindings PRIVATE src/HashingAlgorithms.cpp)
    #    list(APPEND py_files examples/hashing.py)
    #else()
    #    target_sources(ActsPythonBindings PRIVATE stubs/HashingAlgorithmsStub.cpp)
    #endif()

    # Additional components built as separate libraries

    # Podio
    #if(ACTS_BUILD_EXAMPLES_PODIO)
    #    pybind11_add_module(ActsPythonBindingsPodio src/PodioComponent.cpp)
    #    target_link_libraries(
    #        ActsPythonBindingsPodio
    #        PUBLIC ActsPythonUtilities ActsExamplesIoPodio
    #    )
    #    add_dependencies(ActsPythonBindings ActsPythonBindingsPodio)
    #    install(TARGETS ActsPythonBindingsPodio DESTINATION ${_python_install_dir})
    #    set_target_properties(
    #        ActsPythonBindingsPodio
    #        PROPERTIES INSTALL_RPATH "\$ORIGIN/../../${CMAKE_INSTALL_LIBDIR}"
    #    )
    #    set_target_properties(
    #        ActsPythonBindingsPodio
    #        PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_python_dir}/acts
    #    )
    #    list(APPEND py_files examples/podio.py)
    #endif()

    # EDM4hep
    #if(ACTS_BUILD_EXAMPLES_EDM4HEP)
    #    pybind11_add_module(ActsPythonBindingsEDM4hep src/EDM4hepComponent.cpp)
    #    target_link_libraries(
    #        ActsPythonBindingsEDM4hep
    #        PUBLIC ActsPythonUtilities ActsExamplesIoEDM4hep
    #    )
    #    add_dependencies(ActsPythonBindings ActsPythonBindingsEDM4hep)
    #    install(
    #        TARGETS ActsPythonBindingsEDM4hep DESTINATION ${_python_install_dir}
    #    )
    #    set_target_properties(
    #        ActsPythonBindingsEDM4hep
    #        PROPERTIES INSTALL_RPATH "\$ORIGIN/../../${CMAKE_INSTALL_LIBDIR}"
    #    )
    #    set_target_properties(
    #        ActsPythonBindingsEDM4hep
    #        PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_python_dir}/acts
    #    )
    #    list(APPEND py_files examples/edm4hep.py)
    #endif()

    # Geant4
    #if(ACTS_BUILD_EXAMPLES_GEANT4)
    #    pybind11_add_module(ActsPythonBindingsExamplesGeant4 src/Geant4Component.cpp)
    #    target_link_libraries(
    #        ActsPythonBindingsExamplesGeant4
    #        PUBLIC
    #        ActsPythonUtilities
    #        ActsExamplesGeant4
    #        ActsExamplesDetectorGeant4
    #        ActsExamplesMuonSpectrometerMockupDetector
    #    )
    #    add_dependencies(ActsPythonBindings ActsPythonBindingsExamplesGeant4)

    #    install(TARGETS ActsPythonBindingsExamplesGeant4 DESTINATION ${_python_install_dir})
    #    set_target_properties(
    #        ActsPythonBindingsExamplesGeant4
    #        PROPERTIES INSTALL_RPATH "\$ORIGIN/../../${CMAKE_INSTALL_LIBDIR}"
    #    )
    #    set_target_properties(
    #        ActsPythonBindingsExamplesGeant4
    #        PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_python_dir}/acts
    #    )
    #    list(APPEND py_files examples/geant4/__init__.py)
    #endif()

    # DD4hep
    #if(ACTS_BUILD_PLUGIN_DD4HEP AND ACTS_BUILD_EXAMPLES_DD4HEP)
    #    pybind11_add_module(ActsPythonBindingsExamplesDD4hep src/DD4hepComponent.cpp)
    #    target_link_libraries(
    #        ActsPythonBindingsExamplesDD4hep
    #        PUBLIC ActsPythonUtilities ActsExamplesDetectorDD4hep
    #    )
    #    add_dependencies(ActsPythonBindings ActsPythonBindingsExamplesDD4hep)

    #    install(TARGETS ActsPythonBindingsExamplesDD4hep DESTINATION ${_python_install_dir})
    #    set_target_properties(
    #        ActsPythonBindingsExamplesDD4hep
    #        PROPERTIES INSTALL_RPATH "\$ORIGIN/../../${CMAKE_INSTALL_LIBDIR}"
    #    )
    #    set_target_properties(
    #        ActsPythonBindingsExamplesDD4hep
    #        PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_python_dir}/acts
    #    )
    #    list(APPEND py_files examples/dd4hep.py)
    #endif()
else()
    target_sources(
        ActsPythonBindings
        PRIVATE
        src/ExamplesModuleEntryStub.cpp)
endif()