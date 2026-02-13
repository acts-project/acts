acts_add_library(ExamplesDispatchAlgorithms src/PatternDispatchAlgorithm.cpp)

target_include_directories(
    ActsExamplesDispatchAlgorithms
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
    ActsExamplesDispatchAlgorithms
    PUBLIC Acts::ExamplesFramework
)
acts_compile_headers(ExamplesDispatchAlgorithms GLOB "include/**/*.hpp")
