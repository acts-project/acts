if(ACTS_BUILD_PLUGIN_GEOMODEL OR ACTS_BUILD_EXAMPLES_GEANT4)
    acts_add_library(ExamplesMuonSpectrometerMockupDetector SHARED)

    target_link_libraries(
        ActsExamplesMuonSpectrometerMockupDetector
        PUBLIC Acts::ExamplesFramework
    )
    target_include_directories(
        ActsExamplesMuonSpectrometerMockupDetector
        PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    )
    if(ACTS_BUILD_PLUGIN_GEOMODEL)
        target_sources(
            ActsExamplesMuonSpectrometerMockupDetector
            PUBLIC src/GeoMuonMockupExperiment.cpp
        )
        target_link_libraries(
            ActsExamplesMuonSpectrometerMockupDetector
            PUBLIC ActsPluginGeoModel GeoModelIO::GeoModelWrite
        )
    endif()
    if(ACTS_BUILD_EXAMPLES_GEANT4)
        target_sources(
            ActsExamplesMuonSpectrometerMockupDetector
            PUBLIC src/MockupSectorBuilder.cpp
        )
        target_link_libraries(
            ActsExamplesMuonSpectrometerMockupDetector
            PUBLIC
                Acts::ExamplesGeant4
                Acts::ExamplesDetectorGeant4
                Acts::PluginGeant4
        )
    endif()
endif()
