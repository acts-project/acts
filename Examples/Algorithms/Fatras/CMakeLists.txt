acts_add_library(ExamplesFatras src/FatrasSimulation.cpp)
target_include_directories(
    ActsExamplesFatras
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesFatras
    PUBLIC Acts::Fatras Acts::ExamplesFramework Acts::ExamplesMagneticField
)
