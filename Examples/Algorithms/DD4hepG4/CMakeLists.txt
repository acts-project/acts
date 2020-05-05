add_library(
  ActsExamplesDD4hepG4 SHARED
  src/DD4hepToG4Svc.cpp
  src/GeoConstruction.cpp)
target_include_directories(
  ActsExamplesDD4hepG4
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesDD4hepG4
  PUBLIC
    ActsExamplesFramework ActsExamplesDetectorDD4hep
    ${Geant4_LIBRARIES})

if(${DD4hep_VERSION} VERSION_LESS 1.11)
  target_include_directories(
    ActsExamplesDD4hepG4
    PUBLIC ${DD4hep_INCLUDE_DIRS})
  target_link_libraries(
    ActsExamplesDD4hepG4
    PUBLIC ${DD4hep_DDCORE_LIBRARY} ${DD4hep_DDG4_LIBRARY})
else()
  target_link_libraries(
    ActsExamplesDD4hepG4
    PUBLIC DD4hep::DDCore DD4hep::DDG4)
endif()

install(
  TARGETS ActsExamplesDD4hepG4
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
