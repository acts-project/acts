acts_add_library(
    PluginMille
    src/Helpers.cpp
    src/ActsToMille.cpp
    ACTS_INCLUDE_FOLDER include/ActsPlugins
)

target_include_directories(
    ActsPluginMille
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(ActsPluginMille PUBLIC ActsAlignment Mille2)

acts_compile_headers(PluginMille GLOB include/**/*.hpp)
