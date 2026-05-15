acts_add_library(
    FatrasGeant4
    SHARED
    src/Geant4Decay.cpp
    src/PDGtoG4Converter.cpp
    src/DummyDetectorConstruction.cpp
    ACTS_INCLUDE_FOLDER include/ActsFatras
)
target_compile_definitions(ActsFatrasGeant4 PUBLIC ${Geant4_DEFINITIONS})
target_include_directories(
    ActsFatrasGeant4
    SYSTEM
    PUBLIC ${Geant4_INCLUDE_DIRS}
)
target_include_directories(
    ActsFatrasGeant4
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(ActsFatrasGeant4 PUBLIC Acts::Fatras ${Geant4_LIBRARIES})
