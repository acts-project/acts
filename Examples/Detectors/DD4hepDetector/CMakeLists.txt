acts_add_library(
    ExamplesDetectorDD4hep
    src/DD4hepDetector.cpp
    src/AlignedDD4hepDetectorElement.cpp
    src/OpenDataDetector.cpp
    ACTS_INCLUDE_FOLDER include/ActsExamples
)

target_include_directories(
    ActsExamplesDetectorDD4hep
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(
    ActsExamplesDetectorDD4hep
    PUBLIC Acts::PluginDD4hep Acts::ExamplesDetectorsCommon
)

if(ACTS_BUILD_EXAMPLES_GEANT4)
    if(${DD4hep_VERSION} VERSION_LESS 1.11)
        target_link_libraries(
            ActsExamplesDetectorDD4hep
            PUBLIC Acts::ExamplesGeant4 ${DD4hep_DDG4_LIBRARY}
        )
    else()
        target_link_libraries(
            ActsExamplesDetectorDD4hep
            PUBLIC Acts::ExamplesGeant4 DD4hep::DDG4
        )
    endif()

    target_sources(
        ActsExamplesDetectorDD4hep
        PRIVATE src/DDG4DetectorConstruction.cpp src/DD4hepDetectorGeant4.cpp
    )
else()
    target_sources(
        ActsExamplesDetectorDD4hep
        PRIVATE src/DD4hepDetectorGeant4Stub.cpp
    )
endif()

if(${DD4hep_VERSION} VERSION_LESS 1.11)
    target_include_directories(
        ActsExamplesDetectorDD4hep
        PUBLIC ${DD4hep_INCLUDE_DIRS}
    )
    target_link_libraries(
        ActsExamplesDetectorDD4hep
        PUBLIC ${DD4hep_DDCORE_LIBRARY}
    )
else()
    target_link_libraries(ActsExamplesDetectorDD4hep PUBLIC DD4hep::DDCore)
endif()

# not sure why this needs to be set, but dd4hep fails to compile otherwise
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR})
set_target_properties(
    ActsExamplesDetectorDD4hep
    PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
