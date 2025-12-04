acts_add_library(
    ExamplesDetectorsCommon
    SHARED
    src/Detector.cpp
    src/AlignmentDecorator.cpp
    src/StructureSelector.cpp
    ACTS_INCLUDE_FOLDER include/ActsExamples
)
target_include_directories(
    ActsExamplesDetectorsCommon
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesDetectorsCommon
    PUBLIC Acts::ExamplesFramework
)

acts_compile_headers(ExamplesDetectorsCommon GLOB "include/**/*.hpp")
