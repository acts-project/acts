acts_add_library(
    ExamplesMaterialMapping
    src/CoreMaterialMapping.cpp
    src/MaterialValidation.cpp
    src/MaterialMapping.cpp
)

target_include_directories(
    ActsExamplesMaterialMapping
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesMaterialMapping
    PUBLIC Acts::ExamplesFramework
)

acts_compile_headers(ExamplesMaterialMapping GLOB "include/**/*.hpp")
