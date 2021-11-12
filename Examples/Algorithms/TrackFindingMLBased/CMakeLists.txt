add_library(
  ActsExamplesTrackFindingMLBased SHARED
  src/TrackFindingMLBasedAlgorithm.cpp
)


target_include_directories(
  ActsExamplesTrackFindingMLBased
  PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
    

target_link_libraries(
  ActsExamplesTrackFindingMLBased
  PUBLIC
    ActsExamplesFramework
)

install(
  TARGETS ActsExamplesTrackFindingMLBased
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})