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
    PUBLIC Acts::Core
    PRIVATE Arrow::arrow_static Parquet::parquet_static
)

# Propagate arrow/parquet headers to direct consumers. The static archives
# themselves stay PRIVATE so they only enter libActsPluginArrow once. This
# library is the sole "arrow island" in ACTS: any other .cpp that uses arrow
# directly must be donated in as an OBJECT library (see
# Examples/Io/{Parquet,Arrow}/CMakeLists.txt).
target_include_directories(
    ActsPluginArrow
    SYSTEM
    PUBLIC
        $<TARGET_PROPERTY:Arrow::arrow_static,INTERFACE_INCLUDE_DIRECTORIES>
        $<TARGET_PROPERTY:Parquet::parquet_static,INTERFACE_INCLUDE_DIRECTORIES>
        >>>>>>>
        conflict
        1
        of
        1
        ends
)

# Hide symbols coming from the arrow/parquet static archives so they can't
# collide with pyarrow's bundled libarrow when both are loaded into the same
# Python process. Our own ACTS code keeps default visibility (ACTS does not
# use an ACTS_EXPORT macro; hiding its symbols would break downstream
# consumers). The load_hidden / --exclude-libs,ALL linker flags strip only
# the archive-supplied symbols from the dynamic symbol table.
if(APPLE)
    target_link_options(
        ActsPluginArrow
        PRIVATE
            "LINKER:-load_hidden,$<TARGET_FILE:Arrow::arrow_static>"
            "LINKER:-load_hidden,$<TARGET_FILE:Parquet::parquet_static>"
    )
else()
    target_link_options(ActsPluginArrow PRIVATE "LINKER:--exclude-libs,ALL")
endif()

acts_compile_headers(PluginArrow GLOB include/**/*.hpp)
