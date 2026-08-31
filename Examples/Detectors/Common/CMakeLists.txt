acts_add_library(
    ExamplesDetectorsCommon
    src/Detector.cpp
    src/AlignmentDecorator.cpp
    src/StructureSelector.cpp
    ACTS_INCLUDE_FOLDER include/ActsExamples
)
target_include_directories(
    ActsExamplesDetectorsCommon
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(
    ActsExamplesDetectorsCommon
    PUBLIC Acts::ExamplesFramework
)
