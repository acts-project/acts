acts_add_library(
    ExamplesGeneratorsPythia8
    ActsExamples/Generators/Pythia8ProcessGenerator.cpp
)
target_include_directories(
    ActsExamplesGeneratorsPythia8
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
target_link_libraries(
    ActsExamplesGeneratorsPythia8
    PUBLIC Acts::ExamplesFramework Acts::ExamplesGenerators
    PRIVATE Pythia8::Pythia8
)
