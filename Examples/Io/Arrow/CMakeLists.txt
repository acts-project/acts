add_library(
    ActsExamplesIoArrow_obj
    OBJECT
    src/ArrowCaloHitOutputConverter.cpp
    src/ArrowParticleOutputConverter.cpp
    src/ArrowSimHitOutputConverter.cpp
    src/ArrowTrackOutputConverter.cpp
    src/ColliderMLRelease1InputConverter.cpp
)

target_include_directories(
    ActsExamplesIoArrow_obj
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        # Genex (not target_link_libraries) for Plugin and sibling OBJECT
        # donor: a link-level dep on ActsPluginArrow would create a cycle
        # since it absorbs this OBJECT lib's .o files.
        $<TARGET_PROPERTY:Acts::ExamplesIoParquet,INTERFACE_INCLUDE_DIRECTORIES>
        $<TARGET_PROPERTY:Acts::PluginArrow,INTERFACE_INCLUDE_DIRECTORIES>
        $<TARGET_PROPERTY:Acts::ExamplesDigitization,INTERFACE_INCLUDE_DIRECTORIES>
        $<TARGET_PROPERTY:Acts::ExamplesIoCsv,INTERFACE_INCLUDE_DIRECTORIES>
)

target_link_libraries(
    ActsExamplesIoArrow_obj
    PRIVATE Acts::ExamplesFramework Acts::ArrowLinkage Acts::ParquetLinkage
)

set_target_properties(
    ActsExamplesIoArrow_obj
    PROPERTIES POSITION_INDEPENDENT_CODE ON
)
if(ACTS_ARROW_ISOLATED)
    set_target_properties(
        ActsExamplesIoArrow_obj
        PROPERTIES CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN YES
    )
endif()

target_sources(
    ActsPluginArrow
    PRIVATE $<TARGET_OBJECTS:ActsExamplesIoArrow_obj>
)
# ColliderMLRelease1InputConverter needs ExamplesDigitization at link time.
target_link_libraries(ActsPluginArrow PRIVATE Acts::ExamplesDigitization)
# ExamplesDigitization is PRIVATE-linked by ActsPluginArrow; add it to the
# export set so the install(EXPORT ActsPluginArrowTargets) resolves.
install(TARGETS ActsExamplesDigitization EXPORT ActsPluginArrowTargets)

acts_add_library(ExamplesIoArrow INTERFACE)
target_include_directories(
    ActsExamplesIoArrow
    INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesIoArrow
    INTERFACE
        Acts::ExamplesIoParquet
        Acts::PluginArrow
        Acts::ExamplesFramework
        Acts::ExamplesDigitization
)

acts_compile_headers(ExamplesIoArrow GLOB include/**/*.hpp)
