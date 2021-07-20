add_executable(
  ActsExampleHepMC3Exe
  HepMC3Example.cpp)
target_link_libraries(
  ActsExampleHepMC3Exe
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
    ActsExamplesFramework ActsExamplesIoHepMC3 ActsExamplesHepMC3 ActsExamplesIoNuclearInteractions
    Boost::program_options)
    
install(
  TARGETS ActsExampleHepMC3Exe ActsExampleHepMC3Histogram
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  
