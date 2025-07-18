acts_add_library(
    ExamplesDetectorTGeo
    SHARED
    src/AlignedTGeoDetectorElement.cpp
    src/TGeoDetector.cpp
    src/TGeoITkModuleSplitter.cpp
)

target_include_directories(
    ActsExamplesDetectorTGeo
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
    ActsExamplesDetectorTGeo
    PUBLIC
        Acts::Core
        Acts::PluginRoot
        Acts::PluginJson
        Acts::ExamplesFramework
        Acts::ExamplesDetectorGeneric
        Acts::ExamplesDetectorsCommon
        Acts::ExamplesITkModuleSplitting
)

acts_compile_headers(ExamplesDetectorTGeo GLOB "include/**/*.hpp")

install(
    TARGETS ActsExamplesDetectorTGeo
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY include/ActsExamples DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
