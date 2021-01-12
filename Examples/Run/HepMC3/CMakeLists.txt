add_executable(
  ActsExampleHepMC3
  HepMC3Example.cpp)
target_link_libraries(
  ActsExampleHepMC3
  PRIVATE
    ActsCore ActsFatras 
    ActsExamplesCommon ActsExamplesFramework ActsExamplesIoHepMC3)

add_executable(
  ActsExampleHepMC3Histogram
  HepMC3HistogramExample.cpp)
target_link_libraries(
  ActsExampleHepMC3Histogram
  PRIVATE
    ActsCore ActsExamplesCommon
    ActsExamplesFramework ActsExamplesIoHepMC3 ActsExamplesHepMC ActsExamplesIoNuclearInteractions
    Boost::program_options)
    
install(
  TARGETS ActsExampleHepMC3 ActsExampleHepMC3Histogram
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  
