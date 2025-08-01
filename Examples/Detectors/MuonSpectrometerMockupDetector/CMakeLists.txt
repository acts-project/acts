set(src_files)
set(link_libs)
if(ACTS_BUILD_PLUGIN_GEOMODEL OR ACTS_BUILD_EXAMPLES_GEANT4)
    if(ACTS_BUILD_EXAMPLES_GEANT4)
        set(src_files ${src_files} src/GeoMuonMockupExperiment.cpp)
        set(link_libs
            ${link_libs}
            ActsExamplesGeant4
            ActsExamplesDetectorGeant4
            ActsPluginGeant4
        )
    endif()
    if(ACTS_BUILD_PLUGIN_GEOMODEL)
        set(src_files ${src_files} src/MockupSectorBuilder.cpp)
        set(link_libs ${link_libs} ActsPluginGeoModel GeoModelIO::GeoModelWrite)
    endif()
endif()

add_library(
    ActsExamplesMuonSpectrometerMockupDetector
    SHARED
    src/MuonSpectrometerMockupDetector.cpp
    ${src_files}
)

target_include_directories(
    ActsExamplesMuonSpectrometerMockupDetector
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesMuonSpectrometerMockupDetector
    PUBLIC ActsCore ActsExamplesFramework ${link_libs}
)

unset(link_libs)
unset(src_files)
install(
    TARGETS ActsExamplesMuonSpectrometerMockupDetector
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
