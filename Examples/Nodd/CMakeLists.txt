acts_add_library(
    ExamplesNodd
    src/DD4hepNoddDetector.cpp
    ACTS_INCLUDE_FOLDER include/ActsExamples
)

target_include_directories(
    ActsExamplesNodd
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(ActsExamplesNodd PUBLIC Acts::ExamplesDetectorDD4hep)

if(ACTS_BUILD_EXAMPLES_GEANT4 AND ACTS_BUILD_PYTHON_BINDINGS)
    install(
        PROGRAMS
            Scripts/nodd_material_recording.py
            Scripts/nodd_geometry_validation.py
        DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()
