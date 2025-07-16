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
    target_link_libraries(
        ActsExamplesMuonSpectrometerMockupDetector
        PUBLIC ActsPluginGeoModel GeoModelIO::GeoModelWrite
    )
    target_sources(
        ActsExamplesMuonSpectrometerMockupDetector
        PUBLIC src/GeoMuonMockupExperiment.cpp
    )
endif()

install(
    TARGETS ActsExamplesMuonSpectrometerMockupDetector
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
