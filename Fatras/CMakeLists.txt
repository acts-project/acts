add_library(
  ActsFatras SHARED
  src/Digitization/Channelizer.cpp
  src/Digitization/DigitizationError.cpp
  src/Digitization/PlanarSurfaceMask.cpp
  src/Digitization/PlanarSurfaceDrift.cpp
  src/EventData/Particle.cpp
  src/EventData/ProcessType.cpp
  src/Kernel/SimulatorError.cpp
  src/Physics/StandardPhysicsLists.cpp
  src/Physics/Decay/Decay.cpp
  src/Physics/Decay/PDGtoG4Converter.cpp
  src/Physics/Decay/G4DetectorConstruction.cpp
  src/Utilities/LandauDistribution.cpp
  src/Utilities/ParticleData.cpp)
target_compile_definitions(
  ActsFatras
  PUBLIC ${Geant4_DEFINITIONS})
target_include_directories(
  ActsFatras
  SYSTEM PUBLIC ${Geant4_INCLUDE_DIRS})
target_compile_features(
  ActsFatras
  PUBLIC cxx_std_17)
target_include_directories(
  ActsFatras
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(
  ActsFatras
  PUBLIC ActsCore ${Geant4_LIBRARIES})

install(
  TARGETS ActsFatras
  EXPORT ActsFatrasTargets
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(
  DIRECTORY include/ActsFatras
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
