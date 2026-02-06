acts_add_library(
    ExamplesDetectorTGeo
    src/AlignedTGeoDetectorElement.cpp
    src/TGeoDetector.cpp
    src/TGeoITkModuleSplitter.cpp
    ACTS_INCLUDE_FOLDER include/ActsExamples
)

target_include_directories(
    ActsExamplesDetectorTGeo
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
    ActsExamplesDetectorTGeo
    PUBLIC
        Acts::PluginRoot
        Acts::PluginJson
        Acts::ExamplesFramework
        Acts::ExamplesDetectorGeneric
        Acts::ExamplesDetectorsCommon
        Acts::ExamplesITkModuleSplitting
)

acts_compile_headers(ExamplesDetectorTGeo GLOB "include/**/*.hpp")
