acts_add_library(ExamplesMagneticField src/ScalableBFieldService.cpp)
target_include_directories(
    ActsExamplesMagneticField
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(ActsExamplesMagneticField PUBLIC Acts::ExamplesFramework)

acts_compile_headers(ExamplesMagneticField GLOB "include/**/*.hpp")
