add_library(ActsExamplesDetectorGeoModel SHARED src/dummy.cpp)

target_include_directories(
    ActsExamplesDetectorGeoModel
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
    ActsExamplesDetectorGeoModel
    PUBLIC
        ActsCore
        ActsExamplesFramework
        ActsPluginGeoModel
        ActsExamplesDetectorsCommon
)

if(ACTS_BUILD_EXAMPLES_GEANT4)
    target_link_libraries(ActsExamplesDetectorDD4hep PUBLIC ActsExamplesGeant4)

    target_sources(
        ActsExamplesDetectorDD4hep
        PUBLIC src/GeoModelGeant4DetectorConstruction.cpp
    )
endif()

install(
    TARGETS ActsExamplesDetectorGeoModel
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
