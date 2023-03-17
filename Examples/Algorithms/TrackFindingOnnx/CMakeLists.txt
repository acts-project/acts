add_library(
  ActsExamplesTrackFindingOnnx SHARED
  src/AmbiguityResolutionMLAlgorithm.cpp
)

target_include_directories(
  ActsExamplesTrackFindingOnnx
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
  ActsExamplesTrackFindingOnnx
  PUBLIC
    ActsCore
    ActsExamplesFramework
)

install(
  TARGETS ActsExamplesTrackFindingOnnx
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
