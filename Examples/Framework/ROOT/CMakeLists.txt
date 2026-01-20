acts_add_library(
    ExamplesFrameworkROOT
    src/ScalingCalibrator.cpp
    src/Helpers.cpp
    ACTS_INCLUDE_FOLDER include/ActsExamples
)

target_include_directories(
    ActsExamplesFrameworkROOT
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
    ActsExamplesFrameworkROOT
    PUBLIC Acts::ExamplesFramework ROOT::Core ROOT::Hist
)
