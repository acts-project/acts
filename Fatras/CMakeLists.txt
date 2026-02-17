acts_add_library(
    Fatras
    src/Digitization/Segmentizer.cpp
    src/Digitization/DigitizationError.cpp
    src/Digitization/PlanarSurfaceMask.cpp
    src/Digitization/PlanarSurfaceDrift.cpp
    src/EventData/Particle.cpp
    src/EventData/ParticleOutcome.cpp
    src/EventData/ProcessType.cpp
    src/Kernel/SimulationError.cpp
    src/Physics/BetheHeitler.cpp
    src/Physics/StandardInteractions.cpp
    src/Utilities/LandauDistribution.cpp
    ACTS_INCLUDE_FOLDER include/ActsFatras
)
target_include_directories(
    ActsFatras
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(ActsFatras PUBLIC Acts::Core)

add_component_if(Geant4 FatrasGeant4 ACTS_BUILD_FATRAS_GEANT4)
