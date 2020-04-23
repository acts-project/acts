add_executable(
  ActsExampleHepMC3
  HepMC3Example.cpp)
target_link_libraries(
  ActsExampleHepMC3
  PRIVATE
    ActsCore
    ActsExamplesFramework ActsExamplesIoHepMC3
    Boost::program_options)

install(
  TARGETS ActsExampleHepMC3
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
