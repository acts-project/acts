acts_add_library(
    ExamplesFrameworkML
    src/NeuralCalibrator.cpp
    ACTS_INCLUDE_FOLDER include/ActsExamples
)

target_include_directories(
    ActsExamplesFrameworkML
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
    ActsExamplesFrameworkML
    PUBLIC Acts::ExamplesFramework Acts::PluginOnnx
)
