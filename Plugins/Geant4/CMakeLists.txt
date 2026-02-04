acts_add_library(
    PluginGeant4
    src/Geant4Converters.cpp
    src/Geant4DetectorElement.cpp
    src/Geant4DetectorSurfaceFactory.cpp
    src/Geant4PhysicalVolumeSelectors.cpp
    ACTS_INCLUDE_FOLDER include/ActsPlugins
)
target_include_directories(
    ActsPluginGeant4
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_include_directories(
    ActsPluginGeant4
    SYSTEM
    PRIVATE ${Geant4_INCLUDE_DIRS}
)
target_link_libraries(ActsPluginGeant4 PUBLIC Acts::Core ${Geant4_LIBRARIES})
