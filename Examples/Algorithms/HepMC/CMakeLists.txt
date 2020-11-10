add_library(
  ActsExamplesProcessExtraction SHARED
  src/HepMCProcessExtractor.cpp)
target_include_directories(
  ActsExamplesProcessExtraction
  SYSTEM PUBLIC ${HEPMC3_INCLUDE_DIR})
target_include_directories(
  ActsExamplesProcessExtraction
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesProcessExtraction
  PUBLIC ActsCore ActsExamplesFramework ActsExamplesIoHepMC3)

install(
  TARGETS ActsExamplesProcessExtraction
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
