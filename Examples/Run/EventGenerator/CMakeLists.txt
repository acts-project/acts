add_executable(
  ActsExampleGenParticleGun
  ParticleGun.cpp)
target_link_libraries(
  ActsExampleGenParticleGun
  PRIVATE
    ActsExamplesFramework ActsExamplesGenerators ActsExamplesCommon
    ActsExamplesIoCsv ActsExamplesIoRoot ActsExamplesPrinters
    Boost::program_options)

install(
  TARGETS ActsExampleGenParticleGun
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_subdirectory_if(Pythia8 ACTS_BUILD_EXAMPLES_PYTHIA8)
