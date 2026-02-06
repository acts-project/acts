acts_add_library(
    ExamplesDetectorsCommon
    src/Detector.cpp
    src/AlignmentDecorator.cpp
    src/StructureSelector.cpp
)
target_include_directories(
    ActsExamplesDetectorsCommon
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesDetectorsCommon
    PUBLIC Acts::ExamplesFramework
)
