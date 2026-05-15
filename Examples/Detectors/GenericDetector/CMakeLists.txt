acts_add_library(
    ExamplesDetectorGeneric
    src/AlignedGenericDetector.cpp
    src/GenericDetector.cpp
    src/GenericDetectorElement.cpp
    src/LayerBuilder.cpp
    src/ProtoLayerCreator.cpp
    src/GenericDetectorBuilder.cpp
)

acts_compile_headers(ExamplesDetectorGeneric GLOB "include/**/*.hpp")

target_include_directories(
    ActsExamplesDetectorGeneric
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
    ActsExamplesDetectorGeneric
    PUBLIC Acts::ExamplesDetectorsCommon
)
