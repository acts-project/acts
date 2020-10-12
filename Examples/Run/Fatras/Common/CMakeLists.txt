add_library(
  ActsExamplesFatrasCommon SHARED
  Fatras.cpp
  FatrasDigitization.cpp
  FatrasInput.cpp
  FatrasSimulation.cpp)
target_include_directories(
  ActsExamplesFatrasCommon
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(
  ActsExamplesFatrasCommon
  PUBLIC
    ActsCore ActsFatras
    ActsExamplesFramework
    ActsExamplesCommon
    ActsExamplesGenerators
    ActsExamplesMagneticField ActsExamplesDetectorsCommon
    ActsExamplesFatras ActsExamplesDigitization
    ActsExamplesIoCsv ActsExamplesIoRoot
    Boost::program_options)

install(
  TARGETS ActsExamplesFatrasCommon
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
