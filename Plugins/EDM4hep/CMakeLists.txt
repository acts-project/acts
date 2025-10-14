acts_add_library(
    PluginEDM4hep
    SHARED
    src/EDM4hepUtil.cpp
    ACTS_INCLUDE_FOLDER include/ActsPlugins
)

target_include_directories(
    ActsPluginEDM4hep
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(
    ActsPluginEDM4hep
    PUBLIC Acts::PluginPodio EDM4HEP::edm4hep
)
