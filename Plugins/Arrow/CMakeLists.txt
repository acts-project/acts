acts_add_library(
    PluginArrow
    src/ArrowUtil.cpp
    ACTS_INCLUDE_FOLDER include/ActsPlugins
)

target_include_directories(
    ActsPluginArrow
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(
    ActsPluginArrow
    PUBLIC Acts::Core Arrow::arrow_shared Parquet::parquet_shared
)

acts_compile_headers(PluginArrow GLOB include/**/*.hpp)
