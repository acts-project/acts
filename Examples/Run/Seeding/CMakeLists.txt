add_executable(
    ActsExampleSeedingGeneric
    GenericSeedingExample.cpp)
target_link_libraries(
    ActsExampleSeedingGeneric
    PRIVATE
    ActsExamplesFramework
    ActsExamplesCommon
    ActsExamplesSeeding
    ActsExamplesDetectorGeneric
    ActsExamplesMagneticField
    ActsExamplesIoCsv
    ActsExamplesIoPerformance)

install(
  TARGETS
  ActsExampleSeedingGeneric
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
