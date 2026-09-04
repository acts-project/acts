acts_add_library(ExamplesGenerators ActsExamples/Generators/EventGenerator.cpp)
target_include_directories(
    ActsExamplesGenerators
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
target_link_libraries(
    ActsExamplesGenerators
    PUBLIC Acts::ExamplesIoHepMC3 HepMC3::HepMC3
)
