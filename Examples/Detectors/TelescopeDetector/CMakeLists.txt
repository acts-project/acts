acts_add_library(
    ExamplesDetectorTelescope
    src/TelescopeDetector.cpp
    src/TelescopeDetectorElement.cpp
    src/BuildTelescopeDetector.cpp
)

target_include_directories(
    ActsExamplesDetectorTelescope
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
    ActsExamplesDetectorTelescope
    PUBLIC Acts::ExamplesFramework Acts::ExamplesDetectorsCommon
)

if(ACTS_BUILD_EXAMPLES_GEANT4)
    target_link_libraries(
        ActsExamplesDetectorTelescope
        PUBLIC Acts::ExamplesGeant4
    )

    target_sources(
        ActsExamplesDetectorTelescope
        PUBLIC
            src/TelescopeDetectorGeant4.cpp
            src/TelescopeG4DetectorConstruction.cpp
    )
else()
    target_sources(
        ActsExamplesDetectorTelescope
        PUBLIC src/TelescopeDetectorGeant4Stub.cpp
    )
endif()
