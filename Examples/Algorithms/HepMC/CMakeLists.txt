add_library(
  ActsExamplesHepMC SHARED
  src/HepMCProcessExtractor.cpp)
target_include_directories(
  ActsExamplesHepMC
  SYSTEM PUBLIC ${HEPMC3_INCLUDE_DIR})
target_include_directories(
  ActsExamplesHepMC
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesHepMC
  PUBLIC ActsCore ActsExamplesFramework ActsExamplesIoHepMC3)

install(
  TARGETS ActsExamplesHepMC
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
