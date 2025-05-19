add_library(
    ActsExamplesMuonSpectrometerMockupDetector
    SHARED
    src/MockupSectorBuilder.cpp
)
target_include_directories(
    ActsExamplesMuonSpectrometerMockupDetector
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesMuonSpectrometerMockupDetector
    PUBLIC
        ActsCore
        ActsExamplesFramework
        ActsExamplesGeant4
        ActsExamplesDetectorGeant4
        ActsPluginGeant4
)
if(ACTS_BUILD_PLUGIN_GEOMODEL)
    add_library(
        ActsExamplesGeoModelMuonMockupBuilder
        SHARED
        src/GeoModelMuonMockupBuilder.cpp
    )

    target_include_directories(
        ActsExamplesGeoModelMuonMockupBuilder
        PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    )

    target_link_libraries(
        ActsExamplesGeoModelMuonMockupBuilder
        PUBLIC
            ActsCore
            ActsExamplesFramework
            ActsExamplesGeant4
            ActsExamplesDetectorGeant4
            ActsPluginGeoModel
    )

    install(
        TARGETS ActsExamplesGeoModelMuonMockupBuilder
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif()

install(
    TARGETS ActsExamplesMuonSpectrometerMockupDetector
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
