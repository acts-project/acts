acts_add_library(ExamplesTraccc
  src/TracccChain.cpp
  src/TracccSeqAlg.cpp
)

target_include_directories(
    ActsExamplesTraccc
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
    ActsExamplesTraccc
    PUBLIC
        Acts::ExamplesFramework
        Acts::ExamplesPropagation
        Acts::PluginDetray
        Acts::PluginCovfie
        traccc::io
        traccc::core
        $<$<BOOL:${ACTS_ENABLE_CUDA}>:traccc::cuda>
        $<$<BOOL:${ACTS_ENABLE_CUDA}>:traccc::device_common>
        $<$<BOOL:${ACTS_ENABLE_CUDA}>:vecmem::cuda>
)

acts_compile_headers(ExamplesTraccc GLOB "include/**/*.hpp")
