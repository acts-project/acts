acts_add_library(
    PluginJson
    src/AlgebraJsonConverter.cpp
    src/ExtentJsonConverter.cpp
    src/AxisSpecJsonConverter.cpp
    src/GridJsonConverter.cpp
    src/GeometryIdentifierJsonConverter.cpp
    src/IndexGridNavigationJsonConverter.cpp
    src/JsonMaterialDecorator.cpp
    src/MaterialMapJsonConverter.cpp
    src/MaterialJsonConverter.cpp
    src/ProtoAxisJsonConverter.cpp
    src/SurfaceBoundsJsonConverter.cpp
    src/SurfaceJsonConverter.cpp
    src/SurfaceMaterialJsonConverter.cpp
    src/UtilitiesJsonConverter.cpp
    src/VolumeBoundsJsonConverter.cpp
    src/VolumeJsonConverter.cpp
    src/AmbiguityConfigJsonConverter.cpp
    src/DetrayJsonHelper.cpp
    src/JsonDetectorElement.cpp
    src/JsonSurfacesReader.cpp
    src/DefinitionsJsonConverter.cpp
    src/SeedingConfigJsonConverter.cpp
    src/TrackingGeometryJsonConverter.cpp
    src/TrackingGeometryMaterialJsonConverter.cpp
    src/BetheHeitlerApproxJsonConverter.cpp
    src/detail/JsonIo.cpp
    ACTS_INCLUDE_FOLDER include/ActsPlugins
)

target_include_directories(
    ActsPluginJson
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(
    ActsPluginJson
    PUBLIC Acts::Core nlohmann_json::nlohmann_json
)

target_compile_definitions(
    ActsPluginJson
    PRIVATE JSON_DISABLE_ENUM_SERIALIZATION=1
)

# zstd is optional: without it, reading a compressed file and requesting
# compressed output both fail with an explanatory error. All zstd use is
# confined to the .cpp, so the definition stays private and the headers
# compile identically either way.
if(ACTS_USE_ZSTD STREQUAL "NEVER")
    message(STATUS "Json plugin: zstd support disabled (ACTS_USE_ZSTD=NEVER)")
else()
    message(
        CHECK_START
        "Json plugin: checking for zstd (ACTS_USE_ZSTD=${ACTS_USE_ZSTD})"
    )
    if(ACTS_USE_ZSTD STREQUAL "ALWAYS")
        find_package(ZSTD REQUIRED)
    else()
        find_package(ZSTD QUIET)
    endif()
    if(ZSTD_FOUND)
        message(CHECK_PASS "found, compression enabled")
        target_link_libraries(ActsPluginJson PRIVATE ZSTD::ZSTD)
        target_compile_definitions(
            ActsPluginJson
            PRIVATE ACTS_JSON_ZSTD_SUPPORT=1
        )
    else()
        message(CHECK_PASS "not found, continuing without compression")
    endif()
endif()

acts_compile_headers(PluginJson GLOB include/**/*.hpp)

# Offline/editor validation; the codec has no runtime schema dependency.
install(
    FILES schema/material-map-v1.schema.json
    DESTINATION ${CMAKE_INSTALL_DATADIR}/Acts/schema
)
