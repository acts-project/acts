acts_add_library(
    ExamplesDigitization
    src/DigitizationAlgorithm.cpp
    src/DigitizationConfig.cpp
    src/DigitizationCoordinatesConverter.cpp
    src/MeasurementCreation.cpp
    src/DigitizationConfigurator.cpp
    src/ModuleClusters.cpp
    src/MuonSpacePointDigitizer.cpp
)
target_include_directories(
    ActsExamplesDigitization
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(ActsExamplesDigitization PUBLIC Acts::ExamplesFramework)

acts_compile_headers(ExamplesDigitization GLOB "include/**/*.hpp")
