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

# Pick static or shared arrow/parquet based on ACTS_ARROW_ISOLATED and
# publish the choice as two INTERFACE targets (Acts::ArrowLinkage and
# Acts::ParquetLinkage) that consumers reference without caring which
# underlying target they forward to. This library is the sole arrow-hosting
# .dylib in ACTS: code that uses arrow directly is donated in via OBJECT
# libraries from Examples/Io/{Parquet,Arrow}. The OBJECT-absorption layout
# does not depend on the isolation toggle.
#
# Include dirs are explicitly mirrored onto the INTERFACE targets (not just
# the link libraries) so downstream uses of INTERFACE_INCLUDE_DIRECTORIES
# resolve without relying on transitive extraction.
if(ACTS_ARROW_ISOLATED)
    set(_acts_arrow_impl Arrow::arrow_static)
    set(_acts_parquet_impl Parquet::parquet_static)
    set(_acts_dataset_impl ArrowDataset::arrow_dataset_static)
else()
    set(_acts_arrow_impl Arrow::arrow_shared)
    set(_acts_parquet_impl Parquet::parquet_shared)
    set(_acts_dataset_impl ArrowDataset::arrow_dataset_shared)
endif()
add_library(ActsArrowLinkage INTERFACE)
add_library(Acts::ArrowLinkage ALIAS ActsArrowLinkage)
target_link_libraries(ActsArrowLinkage INTERFACE ${_acts_arrow_impl})
target_include_directories(
    ActsArrowLinkage
    SYSTEM
    INTERFACE
        $<TARGET_PROPERTY:${_acts_arrow_impl},INTERFACE_INCLUDE_DIRECTORIES>
)
add_library(ActsParquetLinkage INTERFACE)
add_library(Acts::ParquetLinkage ALIAS ActsParquetLinkage)
target_link_libraries(ActsParquetLinkage INTERFACE ${_acts_parquet_impl})
target_include_directories(
    ActsParquetLinkage
    SYSTEM
    INTERFACE
        $<TARGET_PROPERTY:${_acts_parquet_impl},INTERFACE_INCLUDE_DIRECTORIES>
)
add_library(ActsDatasetLinkage INTERFACE)
add_library(Acts::DatasetLinkage ALIAS ActsDatasetLinkage)
target_link_libraries(ActsDatasetLinkage INTERFACE ${_acts_dataset_impl})
target_include_directories(
    ActsDatasetLinkage
    SYSTEM
    INTERFACE
        $<TARGET_PROPERTY:${_acts_dataset_impl},INTERFACE_INCLUDE_DIRECTORIES>
)
# ActsPluginArrow PRIVATE-links these aliases, so CMake requires them to be
# in the same export set for install(EXPORT ActsPluginArrowTargets).
install(
    TARGETS ActsArrowLinkage ActsParquetLinkage ActsDatasetLinkage
    EXPORT ActsPluginArrowTargets
)

target_link_libraries(
    ActsPluginArrow
    PUBLIC Acts::Core
    PRIVATE Acts::ArrowLinkage Acts::ParquetLinkage Acts::DatasetLinkage
)

# Propagate arrow/parquet headers to direct consumers. The linkage itself
# stays PRIVATE on ActsPluginArrow so the archives (isolated mode) or
# libarrow/libparquet DT_NEEDED (shared mode) stay confined to this .dylib
# and don't leak into downstream targets' link lines.
target_include_directories(
    ActsPluginArrow
    SYSTEM
    PUBLIC
        $<TARGET_PROPERTY:Acts::ArrowLinkage,INTERFACE_INCLUDE_DIRECTORIES>
        $<TARGET_PROPERTY:Acts::ParquetLinkage,INTERFACE_INCLUDE_DIRECTORIES>
        $<TARGET_PROPERTY:Acts::DatasetLinkage,INTERFACE_INCLUDE_DIRECTORIES>
)

if(ACTS_ARROW_ISOLATED)
    # Isolate arrow's symbols from every other arrow consumer in the same
    # process (e.g. pyarrow's bundled libarrow) so they cannot collide under
    # any dlopen mode, including Linux RTLD_GLOBAL (e.g. with PyROOT). Three
    # layers of defense, in order of authority:
    #
    # 1. Compile-time visibility: -fvisibility=hidden and
    #    -fvisibility-inlines-hidden mean .o files only mark symbols visible
    #    when explicitly annotated (ACTS_ARROW_EXPORT on the public API).
    # 2. Archive hiding: load_hidden (macOS) / --exclude-libs,ALL (Linux)
    #    strip symbols that came from arrow/parquet's static archives
    #    regardless of how they were emitted.
    # 3. Export allowlist: -exported_symbols_list (macOS) / --version-script
    #    (Linux) drops any symbol whose mangled name doesn't match the
    #    *Acts* pattern. Final authority at link time; catches anything
    #    the first two layers missed (e.g. ARROW_EXPORT-marked template
    #    instantiations in our own TUs).
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
                "LINKER:-load_hidden,$<TARGET_FILE:ArrowDataset::arrow_dataset_static>"
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
endif()

acts_compile_headers(PluginArrow GLOB include/**/*.hpp)
