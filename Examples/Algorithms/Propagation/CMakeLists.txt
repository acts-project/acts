acts_add_library(ExamplesPropagation src/PropagationAlgorithm.cpp)

target_include_directories(
    ActsExamplesPropagation
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(ActsExamplesPropagation PUBLIC Acts::ExamplesFramework)

acts_compile_headers(ExamplesPropagation GLOB "include/**/*.hpp")
