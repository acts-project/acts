acts_add_library(
    PluginHashing
    INTERFACE
    ACTS_INCLUDE_FOLDER include/ActsPlugins
)
target_include_directories(
    ActsPluginHashing
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(ActsPluginHashing INTERFACE Annoy)
