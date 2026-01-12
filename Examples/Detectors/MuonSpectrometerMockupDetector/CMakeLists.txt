if(ACTS_BUILD_PLUGIN_GEOMODEL)
    acts_add_library(ExamplesMuonSpectrometerMockupDetector)

    target_link_libraries(
        ActsExamplesMuonSpectrometerMockupDetector
        PUBLIC Acts::ExamplesFramework
    )
    target_include_directories(
        ActsExamplesMuonSpectrometerMockupDetector
        PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    )
    target_sources(
        ActsExamplesMuonSpectrometerMockupDetector
        PUBLIC src/GeoMuonMockupExperiment.cpp
    )
    target_link_libraries(
        ActsExamplesMuonSpectrometerMockupDetector
        PUBLIC ActsPluginGeoModel GeoModelIO::GeoModelWrite
    )
endif()
