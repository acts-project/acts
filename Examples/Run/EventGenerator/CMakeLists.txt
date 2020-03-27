add_executable(
  ActsGenParticleGun
  ParticleGun.cpp)
target_link_libraries(
  ActsGenParticleGun
  PRIVATE
    ActsExamplesFramework ActsExamplesCommon
    ActsExamplesGenerators ActsExamplesIoRoot ActsExamplesIoCsv
    ActsExamplesPrinters
    Boost::program_options)

install(
  TARGETS ActsGenParticleGun
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_subdirectory_if(Pythia8 ACTS_BUILD_EXAMPLES_PYTHIA8)
