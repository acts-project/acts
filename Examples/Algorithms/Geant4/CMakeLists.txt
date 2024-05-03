set(ACTS_EXAMPLES_G4SOURCES
    src/GdmlDetectorConstruction.cpp
    src/TelescopeG4DetectorConstruction.cpp
    src/Geant4Simulation.cpp
    src/MagneticFieldWrapper.cpp
    src/MaterialPhysicsList.cpp
    src/MaterialSteppingAction.cpp
    src/ParticleTrackingAction.cpp
    src/RegionCreator.cpp
    src/SensitiveSurfaceMapper.cpp
    src/SensitiveSteppingAction.cpp
    src/SimParticleTranslation.cpp
    src/ParticleKillAction.cpp
    src/PhysicsListFactory.cpp
    src/Geant4Manager.cpp)

if (ACTS_BUILD_EXAMPLES_DD4HEP)
    list(APPEND ACTS_EXAMPLES_G4SOURCES src/DDG4DetectorConstruction.cpp)
endif()

add_library(
  ActsExamplesGeant4 SHARED ${ACTS_EXAMPLES_G4SOURCES})
target_compile_definitions(
  ActsExamplesGeant4
  PUBLIC ${Geant4_DEFINITIONS})
target_include_directories(
  ActsExamplesGeant4
  SYSTEM PUBLIC ${Geant4_INCLUDE_DIRS})
target_include_directories(
  ActsExamplesGeant4
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesGeant4
  PUBLIC ActsCore ActsExamplesFramework ActsExamplesDetectorTelescope Boost::headers ${Geant4_LIBRARIES})

if (ACTS_BUILD_EXAMPLES_DD4HEP)
  if(${DD4hep_VERSION} VERSION_LESS 1.11)
  target_include_directories(
    ActsExamplesGeant4
    PRIVATE ${DD4hep_INCLUDE_DIRS})
  target_link_libraries(
    ActsExamplesGeant4
    PRIVATE ${DD4hep_DDCORE_LIBRARY} ${DD4hep_DDG4_LIBRARY})
  else()
  target_link_libraries(
    ActsExamplesGeant4
    PUBLIC ActsExamplesDetectorDD4hep DD4hep::DDCore DD4hep::DDG4)
  endif()
endif()

install(
  TARGETS ActsExamplesGeant4
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
