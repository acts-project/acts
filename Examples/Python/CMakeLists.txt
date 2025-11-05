set(_python_dir "${CMAKE_BINARY_DIR}/python")
file(MAKE_DIRECTORY ${_python_dir})
file(MAKE_DIRECTORY ${_python_dir}/acts)

set(_python_install_dir "python/acts")

target_sources(
    ActsPythonBindings
    PRIVATE
        src/ModuleEntry.cpp
        src/Framework.cpp
        src/Alignment.cpp
        src/Detector.cpp
        src/Digitization.cpp
        src/MaterialMapping.cpp
        src/ExampleAlgorithms.cpp
        src/Output.cpp
        src/Input.cpp
        src/Propagation.cpp
        src/Generators.cpp
        src/Obj.cpp
        src/TruthTracking.cpp
        src/TrackFitting.cpp
        src/TrackFinding.cpp
        src/Vertexing.cpp
        src/AmbiguityResolution.cpp
        src/Utilities.cpp
)

install(TARGETS ActsPythonBindings DESTINATION ${_python_install_dir})

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
    PUBLIC
        Acts::PythonUtilities
        Acts::ExamplesFramework
        Acts::ExamplesGeometry
        Acts::ExamplesDetectorGeneric
        Acts::ExamplesIoCsv
        Acts::ExamplesIoObj
        Acts::ExamplesGenerators
        Acts::ExamplesPrinters
        Acts::ExamplesTrackFinding
        Acts::ExamplesTrackFitting
        Acts::ExamplesVertexing
        Acts::ExamplesFatras
        Acts::ExamplesDetectorTelescope
        Acts::ExamplesUtilities
        Acts::ExamplesAmbiguityResolution
        Acts::ExamplesTruthTracking
        Acts::ExamplesDigitization
        Acts::ExamplesPropagation
        Acts::ExamplesMaterialMapping
)

set(py_files
    examples/__init__.py
    examples/simulation.py
    examples/reconstruction.py
    examples/itk.py
    examples/odd.py
)

#target_link_libraries(ActsPythonBindings PUBLIC Acts::ExamplesIoHepMC3)
#target_sources(ActsPythonBindings PRIVATE src/HepMC3.cpp)
#list(APPEND py_files examples/hepmc3/__init__.py)
#list(APPEND py_files examples/hepmc3/__main__.py)

# Symlink hepmc3_normalize.py to build/bin for development
add_custom_target(
    hepmc_normalize_symlink
    ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/bin
    COMMAND
        ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/python/acts/examples/hepmc_normalize.py
        ${CMAKE_BINARY_DIR}/bin/hepmc_normalize.py
    COMMENT "Creating symlink for hepmc_normalize.py"
    VERBATIM
)
add_dependencies(hepmc_normalize_symlink ActsExamplesIoHepMC3)
install(
    PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/python/acts/examples/hepmc_normalize.py
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    RENAME hepmc_normalize.py
)

foreach(f ${py_files})
    set(_target ${_python_dir}/acts/${f})
    get_filename_component(_dir ${_target} DIRECTORY)
    file(MAKE_DIRECTORY ${_dir})

    file(
        CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/python/acts/${f} ${_target}
        SYMBOLIC
    )

    get_filename_component(_rel ${f} DIRECTORY)

    install(
        FILES ${CMAKE_CURRENT_SOURCE_DIR}/python/acts/${f}
        DESTINATION ${_python_install_dir}/${_rel}
    )
endforeach()
