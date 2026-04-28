acts_add_library(ExamplesTraccc src/TracccChain.cpp src/TracccSeqAlg.cpp)

target_include_directories(
    ActsExamplesTraccc
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
    ActsExamplesTraccc
    PUBLIC
        Acts::ExamplesFramework
        Acts::ExamplesPropagation
        Acts::PluginCovfie
        Acts::PluginDetray
        traccc::io
        traccc::core
        traccc::cuda
        traccc::device_common
        vecmem::cuda
)