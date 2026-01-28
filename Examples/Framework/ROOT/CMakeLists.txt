acts_add_library(
    ExamplesFrameworkRoot
    src/ScalingCalibrator.cpp
    ACTS_INCLUDE_FOLDER include/ActsExamples
)

target_include_directories(
    ActsExamplesFrameworkRoot
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
    ActsExamplesFrameworkRoot
    PUBLIC Acts::ExamplesFramework ROOT::Core ROOT::Hist
)
