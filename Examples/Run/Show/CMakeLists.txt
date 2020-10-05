add_executable(
  ActsExampleShowFatrasGeneric
  ShowFatrasGeneric.cpp)
target_link_libraries(
  ActsExampleShowFatrasGeneric
  PRIVATE
    ActsExamplesFramework ActsExamplesCommon ActsExamplesIoCsv
    ActsExamplesPrinters ActsExamplesDetectorGeneric
    Boost::program_options)
add_executable(
  ActsExampleShowParticles
  ShowParticles.cpp)
target_link_libraries(
  ActsExampleShowParticles
  PRIVATE
    ActsExamplesFramework ActsExamplesCommon ActsExamplesIoCsv
    ActsExamplesPrinters Boost::program_options)

install(
  TARGETS
    ActsExampleShowFatrasGeneric
    ActsExampleShowParticles
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
