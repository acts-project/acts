acts_add_library(
    ExamplesDetectorGeoModel
    src/GeoModelDetector.cpp
    src/AlignedGeoModelDetectorElement.cpp
    src/GeoModelMuonMockupBuilder.cpp
)

target_include_directories(
    ActsExamplesDetectorGeoModel
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
    ActsExamplesDetectorGeoModel
    PUBLIC
        Acts::ExamplesFramework
        Acts::PluginGeoModel
        Acts::ExamplesDetectorsCommon
)

if(ACTS_BUILD_EXAMPLES_GEANT4)
    find_library(GeoModel2G4_LIBRARY GeoModel2G4 REQUIRED)
    # GeoModel >= 6.30 inlines VolumeBuilder's constructor, which reaches
    # Geo2G4MaterialFactory in GeoMaterial2G4 -- a library GeoModel2G4 only
    # links PRIVATE, so it never propagates to us. Not REQUIRED: on 6.27 the
    # symbol still comes from GeoModel2G4 itself.
    find_library(GeoMaterial2G4_LIBRARY GeoMaterial2G4)

    target_link_libraries(
        ActsExamplesDetectorGeoModel
        PUBLIC Acts::ExamplesGeant4 ${GeoModel2G4_LIBRARY}
    )

    if(GeoMaterial2G4_LIBRARY)
        target_link_libraries(
            ActsExamplesDetectorGeoModel
            PUBLIC ${GeoMaterial2G4_LIBRARY}
        )
    endif()

    target_sources(
        ActsExamplesDetectorGeoModel
        PUBLIC
            src/GeoModelGeant4DetectorConstruction.cpp
            src/GeoModelDetectorGeant4.cpp
    )
else()
    target_sources(
        ActsExamplesDetectorGeoModel
        PUBLIC src/GeoModelDetectorGeant4Stub.cpp
    )
endif()
