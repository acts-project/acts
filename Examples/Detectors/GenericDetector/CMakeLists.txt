add_library(
    ActsExamplesDetectorGeneric
    SHARED
    src/GenericDetector.cpp
    src/GenericDetectorElement.cpp
    src/LayerBuilder.cpp
    src/ProtoLayerCreator.cpp
    src/GenericDetectorBuilder.cpp
)

acts_compile_headers(ActsExamplesDetectorGeneric GLOB "include/**/*.hpp")

target_include_directories(
    ActsExamplesDetectorGeneric
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
    ActsExamplesDetectorGeneric
    PUBLIC ActsCore ActsExamplesFramework ActsExamplesDetectorsCommon
)

install(
    TARGETS ActsExamplesDetectorGeneric
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
