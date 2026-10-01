acts_add_library(
    ExamplesIoJson
    src/JsonGeometryList.cpp
    src/JsonMaterialWriter.cpp
    src/JsonSurfacesWriter.cpp
    src/JsonDigitizationConfig.cpp
    ACTS_INCLUDE_FOLDER include/ActsExamples
)
target_include_directories(
    ActsExamplesIoJson
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(
    ActsExamplesIoJson
    PUBLIC
        Acts::PluginJson
        Acts::ExamplesDigitization
        Acts::ExamplesFramework
        Acts::ExamplesMaterialMapping
        Acts::ExamplesTrackFinding
)

acts_compile_headers(ExamplesIoJson GLOB "include/**/*.hpp")
