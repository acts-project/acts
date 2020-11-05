add_library(
  ActsExamplesEventExtraction SHARED
  src/EventExtraction.cpp)
target_include_directories(
  ActsExamplesEventExtraction
  SYSTEM PUBLIC ${HEPMC3_INCLUDE_DIR})
target_include_directories(
  ActsExamplesEventExtraction
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesEventExtraction
  PUBLIC ActsCore ActsExamplesFramework ActsExamplesIoHepMC3 ${HEPMC3_LIBRARIES})

install(
  TARGETS ActsExamplesEventExtraction
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
