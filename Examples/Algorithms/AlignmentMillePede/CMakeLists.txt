acts_add_library(
    ExamplesAlignmentMillePede
    src/MillePedeAlignmentSandbox.cpp
    src/ActsSolverFromMille.cpp
)
target_include_directories(
    ActsExamplesAlignmentMillePede
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesAlignmentMillePede
    PUBLIC
        Acts::Alignment
        Acts::ExamplesFramework
        Acts::ExamplesMagneticField
        ActsPluginMille
)
