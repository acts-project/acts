acts_add_library(
    PluginFastJet
    src/Jets.cpp
    ACTS_INCLUDE_FOLDER include/ActsPlugins
)

target_include_directories(
    ActsPluginFastJet
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(
    ActsPluginFastJet
    PUBLIC Acts::Fatras Acts::Core FastJet::FastJet
)
