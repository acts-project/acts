acts_add_library(
    PluginDetray
    src/DetrayConversionUtils.cpp
    src/DetrayPayloadConverter.cpp
    src/DetrayMaterial.cpp
    src/DetrayNavigation.cpp
    ACTS_INCLUDE_FOLDER include/ActsPlugins
)

# ── Per-(metadata, operation) explicit-instantiation translation units ──────
#
# The heavy detray detector operations (convert, consistency check, JSON
# read/write) are explicitly instantiated for the closed set of metadata types,
# one operation per generated translation unit, so that each TU builds a single
# `detector<metadata>` and stays small and parallel-compilable. Consumers see
# the matching `extern template` declarations and only link these TUs.
#
# Adding a metadata type: add the alias + ACTS_DETRAY_METADATA_FOR_EACH entry in
# include/ActsPlugins/Detray/DetrayMetadata.hpp, then one
# detray_add_metadata_instantiations() call below.
# Adding an operation: add a src/<Op>.cpp.in template and its extern-template
# declaration, then list <Op> here.
set(_detray_instantiation_operations
    DetrayConvert
    DetrayCheckConsistency
    DetrayWriteJson
    DetrayReadDetector
)

function(detray_add_metadata_instantiations META_NAME META_TYPE)
    foreach(op IN LISTS _detray_instantiation_operations)
        set(_generated
            ${CMAKE_CURRENT_BINARY_DIR}/generated/${op}${META_NAME}.cpp
        )
        # @DETRAY_METADATA_TYPE@ is substituted into the template.
        set(DETRAY_METADATA_TYPE ${META_TYPE})
        configure_file(
            ${CMAKE_CURRENT_SOURCE_DIR}/src/${op}.cpp.in
            ${_generated}
            @ONLY
        )
        target_sources(ActsPluginDetray PRIVATE ${_generated})
    endforeach()
endfunction()

detray_add_metadata_instantiations(Odd ActsPlugins::DetrayMetadata::Odd)
detray_add_metadata_instantiations(Default ActsPlugins::DetrayMetadata::Default)

add_dependencies(ActsPluginDetray detray::core covfie::core vecmem::core)

target_include_directories(
    ActsPluginDetray
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(
    ActsPluginDetray
    PUBLIC
        Acts::PluginJson
        detray::core
        detray::core_array
        detray::io
        detray::detectors
        detray::svgtools
        detray::test_utils
        vecmem::core
)

acts_compile_headers(PluginDetray GLOB include/**/*.hpp)
