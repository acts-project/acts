add_executable(
  ActsExampleHepMC3
  HepMC3Example.cpp)
target_link_libraries(
  ActsExampleHepMC3
  PRIVATE
    ActsCore ActsExamplesCommon
    ActsExamplesFramework ActsExamplesIoHepMC3
    Boost::program_options)

install(
  TARGETS ActsExampleHepMC3
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_subdirectory_if(DD4hep ACTS_BUILD_EXAMPLES_DD4HEP AND ACTS_BUILD_EXAMPLES_GEANT4)