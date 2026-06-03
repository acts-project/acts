# The parquet sources use arrow APIs directly. We compile them as an OBJECT
# library and donate the .o files into libActsPluginArrow via
# target_sources(... $<TARGET_OBJECTS:...>). No libActsExamplesIoParquet.so
# is produced; Acts::ExamplesIoParquet becomes an INTERFACE alias that
# forwards to Acts::PluginArrow.
#
# This layout is independent of ACTS_ARROW_ISOLATED — the OBJECT absorption
# is how we keep all arrow-touching code in one .dylib regardless of whether
# arrow is statically or dynamically linked. Only the visibility settings on
# the OBJECT compile and the linkage of arrow itself change with the toggle.
add_library(
    ActsExamplesIoParquet_obj
    OBJECT
    src/ParquetReader.cpp
    src/ParquetWriter.cpp
    src/ArrowInputConverter.cpp
)

target_include_directories(
    ActsExamplesIoParquet_obj
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        # Acts::PluginArrow's includes via genex, NOT target_link_libraries,
        # because a link-level dep would create a cycle: ActsPluginArrow
        # absorbs this OBJECT lib's .o files.
        $<TARGET_PROPERTY:Acts::PluginArrow,INTERFACE_INCLUDE_DIRECTORIES>
)

# Arrow/parquet headers come via the ArrowLinkage/ParquetLinkage aliases
# (defined in Plugins/Arrow/CMakeLists.txt). OBJECT libs don't actually link
# — target_link_libraries here just propagates INTERFACE_INCLUDE_DIRECTORIES
# to the compile step.
target_link_libraries(
    ActsExamplesIoParquet_obj
    PRIVATE
        Acts::ExamplesFramework
        Acts::ArrowLinkage
        Acts::ParquetLinkage
        Acts::DatasetLinkage
)

set_target_properties(
    ActsExamplesIoParquet_obj
    PROPERTIES POSITION_INDEPENDENT_CODE ON
)
if(ACTS_ARROW_ISOLATED)
    set_target_properties(
        ActsExamplesIoParquet_obj
        PROPERTIES CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN YES
    )
endif()

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

# ColliderML input converter — absorbed into libActsPluginArrow like the
# other Parquet TUs. It needs ExamplesDigitization (for createMeasurement and
# smearing types), so we add that as a PRIVATE link to ActsPluginArrow and
# include its headers via genex (not a link-level dep on the OBJECT target).
add_library(
    ActsExamplesIoColliderML_obj
    OBJECT
    src/ColliderMLInputConverter.cpp
)
target_include_directories(
    ActsExamplesIoColliderML_obj
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<TARGET_PROPERTY:Acts::PluginArrow,INTERFACE_INCLUDE_DIRECTORIES>
        $<TARGET_PROPERTY:Acts::ExamplesDigitization,INTERFACE_INCLUDE_DIRECTORIES>
)
target_link_libraries(
    ActsExamplesIoColliderML_obj
    PRIVATE Acts::ExamplesFramework Acts::ArrowLinkage
)
set_target_properties(
    ActsExamplesIoColliderML_obj
    PROPERTIES POSITION_INDEPENDENT_CODE ON
)
if(ACTS_ARROW_ISOLATED)
    set_target_properties(
        ActsExamplesIoColliderML_obj
        PROPERTIES CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN YES
    )
endif()

target_sources(
    ActsPluginArrow
    PRIVATE $<TARGET_OBJECTS:ActsExamplesIoColliderML_obj>
)
target_link_libraries(ActsPluginArrow PRIVATE Acts::ExamplesDigitization)
# ExamplesDigitization is PRIVATE-linked by ActsPluginArrow; add it to the
# export set for the same reason ExamplesFramework is exported above.
install(TARGETS ActsExamplesDigitization EXPORT ActsPluginArrowTargets)

# INTERFACE alias so downstream code using Acts::ExamplesIoColliderML still
# resolves to libActsPluginArrow (no separate .so produced).
acts_add_library(ExamplesIoColliderML INTERFACE)
target_include_directories(
    ActsExamplesIoColliderML
    INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesIoColliderML
    INTERFACE
        Acts::PluginArrow
        Acts::ExamplesFramework
        Acts::ExamplesDigitization
)
