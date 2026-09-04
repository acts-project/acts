acts_add_library(
    ExamplesAlignment
    src/AlignmentAlgorithm.cpp
    src/AlignmentAlgorithmFunction.cpp
)
target_include_directories(
    ActsExamplesAlignment
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesAlignment
    PUBLIC Acts::Alignment Acts::ExamplesFramework Acts::ExamplesMagneticField
)
