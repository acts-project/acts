add_library(
  ActsExamplesFatrasCommon SHARED
  FatrasDigitizationBase.cpp
  FatrasEvgenBase.cpp
  FatrasMain.cpp
  FatrasSimulationBase.cpp)
target_include_directories(
  ActsExamplesFatrasCommon
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(
  ActsExamplesFatrasCommon
  PUBLIC
    ActsCore ActsFatras
    ActsExamplesFramework
    ActsExamplesCommon
    ActsExamplesGenerators ActsExamplesGeneratorsPythia8
    ActsExamplesMagneticField ActsExamplesDetectorsCommon
    ActsExamplesFatras ActsExamplesDigitization
    ActsExamplesIoCsv ActsExamplesIoRoot
    Boost::program_options)

install(
  TARGETS ActsExamplesFatrasCommon
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
