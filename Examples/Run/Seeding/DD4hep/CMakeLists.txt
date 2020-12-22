add_executable(
  ActsExampleSeedingDD4hep
  DD4hepSeedingExample.cpp)

target_link_libraries(
  ActsExampleSeedingDD4hep
  PRIVATE
    ActsExamplesSeedingCommon 
    ActsExamplesDetectorDD4hep)

install(
  TARGETS
    ActsExampleSeedingDD4hep
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
