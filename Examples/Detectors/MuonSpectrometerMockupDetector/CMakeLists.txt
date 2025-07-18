acts_add_library(
    ExamplesMuonSpectrometerMockupDetector
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
        Acts::Core
        Acts::ExamplesFramework
        Acts::ExamplesGeant4
        Acts::ExamplesDetectorGeant4
        Acts::PluginGeant4
)

if(ACTS_BUILD_PLUGIN_GEOMODEL)
    target_link_libraries(
        ActsExamplesMuonSpectrometerMockupDetector
        PUBLIC Acts::PluginGeoModel GeoModelIO::GeoModelWrite
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
