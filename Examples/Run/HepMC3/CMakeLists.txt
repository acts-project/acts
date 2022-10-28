add_library(
  ActsExamplesHepMC3Options SHARED
  Options/src/HepMC3Options.cpp
)
target_include_directories(
  ActsExamplesHepMC3Options PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Options/include>)
target_link_libraries(
  ActsExamplesHepMC3Options
  PRIVATE
    ActsExamplesFramework
    ActsExamplesIoHepMC3
    Boost::program_options)

add_executable(
  ActsExampleHepMC3Exe
  HepMC3Example.cpp)
target_link_libraries(
  ActsExampleHepMC3Exe
  PRIVATE
    ActsCore ActsFatras 
    ActsExamplesCommon ActsExamplesFramework ActsExamplesIoHepMC3 ActsExamplesHepMC3Options)

add_executable(
  ActsExampleHepMC3Histogram
  HepMC3HistogramExample.cpp)
target_link_libraries(
  ActsExampleHepMC3Histogram
  PRIVATE
    ActsCore ActsExamplesCommon
    ActsExamplesFramework ActsExamplesIoHepMC3 ActsExamplesHepMC3 ActsExamplesIoNuclearInteractions ActsExamplesHepMC3Options
    Boost::program_options)
    
install(
  TARGETS ActsExampleHepMC3Exe ActsExampleHepMC3Histogram
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  
