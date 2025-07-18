acts_add_library(
    ExamplesIoPodio
    SHARED
    src/PodioWriter.cpp
    src/PodioReader.cpp
    src/CollectionBaseWriteHandle.cpp
)

target_include_directories(
    ActsExamplesIoPodio
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesIoPodio
    PUBLIC podio::podio Acts::ExamplesFramework Acts::PluginPodio Acts::Core
)

install(TARGETS ActsExamplesIoPodio LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
