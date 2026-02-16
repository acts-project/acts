acts_add_library(
    ExamplesMaterialMapping
    src/CoreMaterialMapping.cpp
    src/MaterialValidation.cpp
    src/MaterialMapping.cpp
    ACTS_INCLUDE_FOLDER include/ActsExamples
)

target_include_directories(
    ActsExamplesMaterialMapping
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(
    ActsExamplesMaterialMapping
    PUBLIC Acts::ExamplesFramework
)

acts_compile_headers(ExamplesMaterialMapping GLOB "include/**/*.hpp")
