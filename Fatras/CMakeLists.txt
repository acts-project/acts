
# if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#     set(_FATRAS_GAMMA_DIST__MACOSX_)
#     message(STATUS "Fatras gamma dist - this looks like MacOS X")
# else()
#     set(_FATRAS_GAMMA_DIST__LINUX_)
#     message(STATUS "Fatras gamma dist - this looks like Linux")
# endif()
# message(STATUS "Did you get any message?")

add_library(
    ActsFatras
    SHARED
    src/Digitization/Segmentizer.cpp
    src/Digitization/DigitizationError.cpp
    src/Digitization/PlanarSurfaceMask.cpp
    src/Digitization/PlanarSurfaceDrift.cpp
    src/EventData/Particle.cpp
    src/EventData/ParticleOutcome.cpp
    src/EventData/ProcessType.cpp
    src/Kernel/SimulationError.cpp
    src/Physics/BetheHeitler.cpp
    src/Physics/PhotonConversion.cpp
    src/Physics/StandardInteractions.cpp
    src/Utilities/LandauDistribution.cpp
)
target_include_directories(
    ActsFatras
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(ActsFatras PUBLIC ActsCore)

install(
    TARGETS ActsFatras
    EXPORT ActsFatrasTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(DIRECTORY include/ActsFatras DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_component_if(Geant4 FatrasGeant4 ACTS_BUILD_FATRAS_GEANT4)
