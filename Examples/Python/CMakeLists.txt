set(_python_dir "${CMAKE_BINARY_DIR}/python")
file(MAKE_DIRECTORY ${_python_dir})
file(MAKE_DIRECTORY ${_python_dir}/acts)

set(_python_install_dir "python/acts")

target_sources(
    ActsPythonBindings
    PRIVATE src/ModuleEntry.cpp src/Propagation.cpp
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
        Acts::ExamplesIoCsv
        Acts::ExamplesIoObj
        Acts::ExamplesGenerators
        Acts::ExamplesPrinters
        Acts::ExamplesTrackFinding
        Acts::ExamplesTrackFitting
        Acts::ExamplesTruthTracking
        Acts::ExamplesVertexing
        Acts::ExamplesFatras
        Acts::ExamplesDigitization
        Acts::ExamplesUtilities
        Acts::ExamplesAmbiguityResolution
        Acts::ExamplesPropagation
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
