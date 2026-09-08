acts_add_library(ExamplesMagneticField
    src/ScalableBFieldService.cpp
    ACTS_INCLUDE_FOLDER include/ActsExamples
)
target_include_directories(
    ActsExamplesMagneticField
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(ActsExamplesMagneticField PUBLIC Acts::ExamplesFramework)

acts_compile_headers(ExamplesMagneticField GLOB "include/**/*.hpp")
