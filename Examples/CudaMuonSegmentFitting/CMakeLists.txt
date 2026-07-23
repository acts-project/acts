acts_add_library(
    ExamplesCudaMuonSegmentFitting
    src/CudaMuonSpacePoint.cpp
    ACTS_INCLUDE_FOLDER include/ActsExamples
)

target_include_directories(
    ActsExamplesCudaMuonSegmentFitting
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
    ActsExamplesCudaMuonSegmentFitting
    PUBLIC Acts::ExamplesFramework CUDA::cudart
)

acts_compile_headers(ExamplesCudaMuonSegmentFitting GLOB "include/**/*.hpp")
