add_library(
    ActsExamplesIoPodio
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
    PUBLIC podio::podio ActsExamplesFramework ActsPluginPodio ActsCore
)

install(TARGETS ActsExamplesIoPodio LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
