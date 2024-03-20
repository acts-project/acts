add_executable(
  ActsExampleShowFatrasGeneric
  ShowFatrasGeneric.cpp)
target_link_libraries(
  ActsExampleShowFatrasGeneric
  PRIVATE
    ActsExamplesFramework ActsExamplesCommon ActsExamplesIoCsv
    ActsExamplesPrinters ActsExamplesDetectorGeneric)
add_executable(
  ActsExampleShowParticles
  ShowParticles.cpp)
target_link_libraries(
  ActsExampleShowParticles
  PRIVATE
    ActsExamplesFramework ActsExamplesCommon ActsExamplesIoCsv
    ActsExamplesPrinters)

install(
  TARGETS
    ActsExampleShowFatrasGeneric
    ActsExampleShowParticles
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
