acts_add_library(
    Fatras
    src/Digitization/Channelizer.cpp
    src/Digitization/Segmentizer.cpp
    src/Digitization/DigitizationError.cpp
    src/Digitization/SurfaceMask.cpp
    src/Digitization/SurfaceDrift.cpp
    src/EventData/Particle.cpp
    src/EventData/SimulationOutcome.cpp
    src/EventData/GenerationProcess.cpp
    src/Kernel/SimulationError.cpp
    src/Physics/BetheHeitler.cpp
    src/Physics/StandardInteractions.cpp
    src/Synthetic/DetectorLayout.cpp
    src/Synthetic/EventConfig.cpp
    src/Synthetic/EventCsvWriter.cpp
    src/Synthetic/EventGenerator.cpp
    src/Synthetic/SeedingTruth.cpp
    src/Synthetic/detail/ParticlePropagator.cpp
    src/Synthetic/detail/PrimaryGenerator.cpp
    src/Synthetic/detail/Propagation.cpp
    src/Synthetic/detail/Sampling.cpp
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

acts_compile_headers(Fatras GLOB include/**/*.hpp)

add_component_if(Geant4 FatrasGeant4 ACTS_BUILD_FATRAS_GEANT4)
add_component_if(Json FatrasJson ACTS_BUILD_PLUGIN_JSON)
