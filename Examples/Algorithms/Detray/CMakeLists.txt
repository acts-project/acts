acts_add_library(ExamplesDetray
  src/ActsToDetrayDetectorAlg.cpp
)

target_include_directories(
  ActsExamplesDetray
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
  ActsExamplesDetray
  PUBLIC
    Acts::ExamplesFramework
    Acts::ExamplesPropagation
    Acts::PluginDetray
    Acts::PluginCovfie
    detray::test_common
)

acts_compile_headers(ExamplesDetray GLOB "include/**/*.hpp")