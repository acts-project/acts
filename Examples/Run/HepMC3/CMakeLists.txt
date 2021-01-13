add_executable(
  ActsExampleHepMC3
  HepMC3Example.cpp)
target_link_libraries(
  ActsExampleHepMC3
  PRIVATE
    ActsCore ActsFatras 
    ActsExamplesCommon ActsExamplesFramework ActsExamplesIoHepMC3)

install(
  TARGETS ActsExampleHepMC3
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
