acts_add_library(
    ExamplesTrackFinding
    src/SpacePointMaker.cpp
    src/TrackFindingAlgorithm.cpp
    src/TrackFindingAlgorithmFunction.cpp
    src/HoughTransformSeeder.cpp
    src/AdaptiveHoughTransformSeeder.cpp
    src/TrackParamsEstimationAlgorithm.cpp
    src/MuonHoughSeeder.cpp
    src/GbtsSeedingAlgorithm.cpp
    src/TrackParamsLookupEstimation.cpp
    src/GridTripletSeedingAlgorithm.cpp
    src/OrthogonalTripletSeedingAlgorithm.cpp
    ACTS_INCLUDE_FOLDER include/ActsExamples
)

target_include_directories(
    ActsExamplesTrackFinding
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(
    ActsExamplesTrackFinding
    PUBLIC Acts::ExamplesMagneticField
)

# If Hashing examples are enabled, add them to the build
if(ACTS_BUILD_EXAMPLES_HASHING)
    target_sources(
        ActsExamplesTrackFinding
        PRIVATE src/HashingPrototypeSeedingAlgorithm.cpp
    )

    target_link_libraries(ActsExamplesTrackFinding PUBLIC Annoy)
endif()

acts_compile_headers(ExamplesTrackFinding GLOB "include/**/*.hpp")
