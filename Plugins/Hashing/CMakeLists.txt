acts_add_library(
    PluginHashing
    SHARED
    src/HashingAlgorithm.cpp
    src/HashingTraining.cpp
    ACTS_INCLUDE_FOLDER include/ActsPlugins
)

target_include_directories(
    ActsPluginHashing
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(ActsPluginHashing PUBLIC Acts::Core Annoy)
