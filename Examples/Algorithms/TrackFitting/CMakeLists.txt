acts_add_library(
    ExamplesTrackFitting
    src/RefittingCalibrator.cpp
    src/TrackFittingAlgorithm.cpp
    src/KalmanFitterFunction.cpp
    src/RefittingAlgorithm.cpp
    src/GsfFitterFunction.cpp
    src/GlobalChiSquareFitterFunction.cpp
)
target_include_directories(
    ActsExamplesTrackFitting
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesTrackFitting
    PUBLIC Acts::ExamplesFramework Acts::ExamplesMagneticField
)

acts_compile_headers(ExamplesTrackFitting GLOB "include/**/*.hpp")
