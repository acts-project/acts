include(${Geant4_USE_FILE})

add_library(
  ActsExamplesDD4hepG4 SHARED
  src/DD4hepToG4Svc.cpp
  src/GeoConstruction.cpp)
target_include_directories(
  ActsExamplesDD4hepG4
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    ${DD4hep_INCLUDE_DIRS})
target_link_libraries(
  ActsExamplesDD4hepG4
  PUBLIC
    ActsExamplesFramework ActsExamplesDetectorDD4hep
    ${DD4hep_LIBRARIES} ${DD4hep_DDG4_LIBRARY} ${Geant4_LIBRARIES})

install(
  TARGETS ActsExamplesDD4hepG4
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
