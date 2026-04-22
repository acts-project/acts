# See Examples/Io/Parquet/CMakeLists.txt for the OBJECT-library rationale.
add_library(ActsExamplesIoArrow_obj OBJECT src/ArrowParticleOutputConverter.cpp)

target_include_directories(
    ActsExamplesIoArrow_obj
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<TARGET_PROPERTY:Acts::ExamplesIoParquet,INTERFACE_INCLUDE_DIRECTORIES>
        $<TARGET_PROPERTY:Acts::PluginArrow,INTERFACE_INCLUDE_DIRECTORIES>
        $<TARGET_PROPERTY:Arrow::arrow_static,INTERFACE_INCLUDE_DIRECTORIES>
        $<TARGET_PROPERTY:Parquet::parquet_static,INTERFACE_INCLUDE_DIRECTORIES>
)

target_link_libraries(ActsExamplesIoArrow_obj PRIVATE Acts::ExamplesFramework)

set_target_properties(
    ActsExamplesIoArrow_obj
    PROPERTIES POSITION_INDEPENDENT_CODE ON
)

target_sources(
    ActsPluginArrow
    PRIVATE $<TARGET_OBJECTS:ActsExamplesIoArrow_obj>
)

acts_add_library(ExamplesIoArrow INTERFACE)
target_include_directories(
    ActsExamplesIoArrow
    INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesIoArrow
    INTERFACE Acts::ExamplesIoParquet Acts::PluginArrow Acts::ExamplesFramework
)

acts_compile_headers(ExamplesIoArrow GLOB include/**/*.hpp)
