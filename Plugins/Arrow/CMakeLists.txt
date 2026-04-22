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

# Isolate arrow's symbols from pyarrow's bundled libarrow (and any other
# arrow consumer in the same process) so they cannot collide under any
# dlopen mode, including Linux RTLD_GLOBAL (e.g. with PyROOT). Three
# layers of defense, in order of authority:
#
# 1. Compile-time visibility: -fvisibility=hidden and
#    -fvisibility-inlines-hidden mean .o files only mark symbols visible
#    when explicitly annotated (ACTS_ARROW_EXPORT on the public API).
# 2. Archive hiding: load_hidden (macOS) / --exclude-libs,ALL (Linux)
#    strip symbols that came from arrow/parquet's static archives
#    regardless of how they were emitted.
# 3. Export allowlist: -exported_symbols_list (macOS) / --version-script
#    (Linux) drop any symbol whose mangled name doesn't match the
#    *Acts* pattern. This is the final authority at link time and catches
#    anything the first two layers missed (e.g. ARROW_EXPORT-marked
#    template instantiations in our own TUs).
set_target_properties(
    ActsPluginArrow
    PROPERTIES CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN YES
)

if(APPLE)
    target_link_options(
        ActsPluginArrow
        PRIVATE
            "LINKER:-load_hidden,$<TARGET_FILE:Arrow::arrow_static>"
            "LINKER:-load_hidden,$<TARGET_FILE:Parquet::parquet_static>"
            "LINKER:-exported_symbols_list,${CMAKE_CURRENT_SOURCE_DIR}/cmake/exported_symbols.txt"
    )
else()
    target_link_options(
        ActsPluginArrow
        PRIVATE
            "LINKER:--exclude-libs,ALL"
            "LINKER:--version-script=${CMAKE_CURRENT_SOURCE_DIR}/cmake/exported_symbols.ld"
    )
endif()

acts_compile_headers(PluginArrow GLOB include/**/*.hpp)
