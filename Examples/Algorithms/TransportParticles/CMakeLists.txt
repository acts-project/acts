add_library(
  ActsExamplesTransportParticles SHARED
  src/TransportParticlesAlgorithm.cpp
)

target_include_directories(
  ActsExamplesTransportParticles
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
  ActsExamplesTransportParticles
  PUBLIC
    ActsCore
    ActsExamplesFramework
    ActsExamplesTrackFitting
    ActsExamplesVertexing
)

install(
  TARGETS ActsExamplesTransportParticles
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
