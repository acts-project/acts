if(ACTS_BUILD_PLUGIN_GEOMODEL AND ACTS_BUILD_EXAMPLES_GEANT4)
    acts_add_library(
        ExamplesMuonSpectrometerMockupDetector
        SHARED
        src/MockupSectorBuilder.cpp
        src/GeoMuonMockupExperiment.cpp
    )

    target_link_libraries(
        ActsExamplesMuonSpectrometerMockupDetector
        PUBLIC
            ActsCore
            ActsExamplesFramework
            ActsExamplesGeant4
            ActsExamplesDetectorGeant4
            ActsPluginGeant4
            ActsPluginGeoModel
            GeoModelIO::GeoModelWrite
    )

    target_include_directories(
        ActsExamplesMuonSpectrometerMockupDetector
        PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    )
    install(
        TARGETS ActsExamplesMuonSpectrometerMockupDetector
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
elseif(ACTS_BUILD_PLUGIN_GEOMODEL)
    acts_add_library(
        ExamplesMuonSpectrometerMockupDetector
        SHARED
        src/GeoMuonMockupExperiment.cpp
    )
    target_link_libraries(
        ActsExamplesMuonSpectrometerMockupDetector
        PUBLIC
            ActsCore
            ActsExamplesFramework
            ActsPluginGeoModel
            GeoModelIO::GeoModelWrite
    )
    target_include_directories(
        ActsExamplesMuonSpectrometerMockupDetector
        PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    )
    install(
        TARGETS ActsExamplesMuonSpectrometerMockupDetector
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
elseif(ACTS_BUILD_PLUGIN_GEOMODEL)
    acts_add_library(
        ExamplesMuonSpectrometerMockupDetector
        SHARED
        src/GeoMuonMockupExperiment.cpp
    )
    target_link_libraries(
        ActsExamplesMuonSpectrometerMockupDetector
        PUBLIC
            ActsCore
            ActsExamplesFramework
            ActsPluginGeant4
            ActsPluginGeoModel
            GeoModelIO::GeoModelWrite
    )
    target_include_directories(
        ActsExamplesMuonSpectrometerMockupDetector
        PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    )
    install(
        TARGETS ActsExamplesMuonSpectrometerMockupDetector
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
elseif(ACTS_BUILD_EXAMPLES_GEANT4)
    acts_add_library(
        ExamplesMuonSpectrometerMockupDetector
        SHARED
        src/MockupSectorBuilder.cpp
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

    target_include_directories(
        ActsExamplesMuonSpectrometerMockupDetector
        PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    )
    install(
        TARGETS ActsExamplesMuonSpectrometerMockupDetector
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif()
