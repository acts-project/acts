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
  src/Utilities/LandauDistribution.cpp
  src/Utilities/ParticleData.cpp)
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
  PUBLIC ActsCore)

install(
  TARGETS ActsFatras
  EXPORT ActsFatrasTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(
  DIRECTORY include/ActsFatras
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  
add_component_if(Geant4 FatrasGeant4 ACTS_BUILD_FATRAS_GEANT4)
