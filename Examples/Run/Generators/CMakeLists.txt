add_executable(
  ActsExampleParticleGun
  ParticleGun.cpp)
target_link_libraries(
  ActsExampleParticleGun
  PRIVATE
    ActsExamplesFramework ActsExamplesGenerators ActsExamplesCommon
    ActsExamplesIoCsv ActsExamplesIoRoot ActsExamplesPrinters)

install(
  TARGETS ActsExampleParticleGun
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_subdirectory_if(Pythia8 ACTS_BUILD_EXAMPLES_PYTHIA8)
