acts_add_library(ExamplesIoArrow src/ArrowParticleOutputConverter.cpp)

acts_compile_headers(ExamplesIoArrow GLOB include/**/*.hpp)

target_include_directories(
    ActsExamplesIoArrow
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesIoArrow
    PUBLIC
        Acts::ExamplesFramework
        Acts::ExamplesIoParquet
        Acts::PluginArrow
        Arrow::arrow_shared
)
