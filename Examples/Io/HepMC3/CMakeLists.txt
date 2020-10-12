add_library(
  ActsExamplesIoHepMC3 SHARED
  src/HepMC3Event.cpp
  src/HepMC3Options.cpp
  src/HepMC3Particle.cpp
  src/HepMC3Reader.cpp
  src/HepMC3Vertex.cpp
  src/HepMC3Writer.cpp
  src/EventRecording.cpp
  src/RunAction.cpp
  src/EventAction.cpp
  src/PrimaryGeneratorAction.cpp
  src/SteppingAction.cpp
  src/HepMC3Options.cpp)
target_include_directories(
  ActsExamplesIoHepMC3
  SYSTEM PUBLIC ${HEPMC3_INCLUDE_DIR})
target_include_directories(
  ActsExamplesIoHepMC3
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesIoHepMC3
  PUBLIC ActsCore ActsExamplesFramework ActsExamplesGeant4 ${HepMC3_LIBRARIES} HepPID)

install(
  TARGETS ActsExamplesIoHepMC3
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
