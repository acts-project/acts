# The parquet sources use arrow APIs directly. To keep arrow's symbols
# confined to a single shared library (libActsPluginArrow, where they are
# statically linked and hidden), we compile these sources as an OBJECT library
# and donate the resulting .o files into libActsPluginArrow via
# target_sources(... $<TARGET_OBJECTS:...>). No libActsExamplesIoParquet.so
# is produced; Acts::ExamplesIoParquet becomes an INTERFACE alias that
# forwards to Acts::PluginArrow.
add_library(
    ActsExamplesIoParquet_obj OBJECT
    src/ParquetReader.cpp
    src/ParquetWriter.cpp
    src/ArrowInputConverter.cpp
)

target_include_directories(
    ActsExamplesIoParquet_obj
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<TARGET_PROPERTY:Acts::PluginArrow,INTERFACE_INCLUDE_DIRECTORIES>
        $<TARGET_PROPERTY:Arrow::arrow_static,INTERFACE_INCLUDE_DIRECTORIES>
        $<TARGET_PROPERTY:Parquet::parquet_static,INTERFACE_INCLUDE_DIRECTORIES>
)

# NOTE: Headers pulled via generator expression on TARGET_PROPERTY, not via
# target_link_libraries(... Acts::PluginArrow). A link-level dep would create
# a cycle since ActsPluginArrow absorbs this OBJECT lib's .o files.
target_link_libraries(
    ActsExamplesIoParquet_obj
    PRIVATE Acts::ExamplesFramework
)

set_target_properties(
    ActsExamplesIoParquet_obj
    PROPERTIES
        POSITION_INDEPENDENT_CODE ON
        CXX_VISIBILITY_PRESET hidden
        VISIBILITY_INLINES_HIDDEN YES
)

target_sources(
    ActsPluginArrow
    PRIVATE $<TARGET_OBJECTS:ActsExamplesIoParquet_obj>
)
# Donating object files doesn't transfer their link deps; propagate the
# non-arrow link requirements of the donated TUs onto libActsPluginArrow.
target_link_libraries(ActsPluginArrow PRIVATE Acts::ExamplesFramework)

# ActsPluginArrow now links ActsExamplesFramework privately, so CMake's
# install(EXPORT ActsPluginArrowTargets) requires ActsExamplesFramework to
# be in the same export set (otherwise: "target ... not in any export set").
# Add it so the export resolves.
install(TARGETS ActsExamplesFramework EXPORT ActsPluginArrowTargets)

# Preserve the Acts::ExamplesIoParquet target for downstream consumers: an
# INTERFACE alias forwarding to libActsPluginArrow, plus the header include
# directory so downstream code can still #include headers from here.
acts_add_library(ExamplesIoParquet INTERFACE)
target_include_directories(
    ActsExamplesIoParquet
    INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesIoParquet
    INTERFACE Acts::PluginArrow Acts::ExamplesFramework
)

acts_compile_headers(ExamplesIoParquet GLOB include/ActsExamples/**/*.hpp)
