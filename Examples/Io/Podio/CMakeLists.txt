acts_add_library(
    ExamplesIoPodio
    src/PodioWriter.cpp
    src/PodioReader.cpp
    src/CollectionBaseWriteHandle.cpp
    src/PodioInputConverter.cpp
    src/PodioMeasurementOutputConverter.cpp
    src/PodioMeasurementInputConverter.cpp
)

target_include_directories(
    ActsExamplesIoPodio
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesIoPodio
    PUBLIC podio::podio Acts::ExamplesFramework Acts::PluginEDM4hep
)
