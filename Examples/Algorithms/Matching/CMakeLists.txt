add_library(
  ActsExamplesMatching SHARED
  src/MatchingAlgorithm.cpp
)

target_include_directories(
  ActsExamplesMatching
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
  ActsExamplesMatching
  PUBLIC
    ActsCore
    ActsExamplesFramework
    ActsExamplesTrackFitting
    ActsExamplesVertexing
)

install(
  TARGETS ActsExamplesMatching
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
