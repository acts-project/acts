# Backend agnostic traccc example code. The backend specific implementations
# are provided by the backend libraries below.
acts_add_library(ExamplesTraccc src/TracccDetector.cpp)

target_include_directories(
    ActsExamplesTraccc
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
    ActsExamplesTraccc
    PUBLIC
        Acts::Core
        Acts::ExamplesFramework
        traccc::core
        vecmem::core
    PRIVATE traccc::io
)

acts_compile_headers(ExamplesTraccc GLOB "include/**/*.hpp")

# CUDA backend
if(ACTS_ENABLE_CUDA)
    acts_add_library(ExamplesTracccCuda Cuda/src/CudaTracccBackend.cpp)

    target_include_directories(
        ActsExamplesTracccCuda
        PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Cuda/include>
    )

    target_link_libraries(
        ActsExamplesTracccCuda
        PUBLIC Acts::ExamplesTraccc
        PRIVATE traccc::cuda vecmem::cuda
    )
endif()
