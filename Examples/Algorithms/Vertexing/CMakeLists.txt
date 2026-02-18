acts_add_library(
    ExamplesVertexing
    src/AdaptiveMultiVertexFinderAlgorithm.cpp
    src/HoughVertexFinderAlgorithm.cpp
    src/IterativeVertexFinderAlgorithm.cpp
    src/VertexFitterAlgorithm.cpp
    src/TruthVertexSeeder.cpp
)
target_include_directories(
    ActsExamplesVertexing
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesVertexing
    PUBLIC Acts::ExamplesFramework
    PRIVATE Acts::ExamplesTruthTracking
)

acts_compile_headers(ExamplesVertexing GLOB "include/**/*.hpp")
