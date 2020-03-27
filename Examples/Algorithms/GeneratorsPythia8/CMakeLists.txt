add_library(
  ActsExamplesGeneratorsPythia8 SHARED
  ACTFW/Generators/Pythia8ProcessGenerator.cpp
  ACTFW/Options/Pythia8Options.cpp)
target_include_directories(
  ActsExamplesGeneratorsPythia8
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(
  ActsExamplesGeneratorsPythia8
  PUBLIC
    ActsCore
    ActsExamplesFramework ActsExamplesGenerators
    Boost::program_options
  PRIVATE
    Pythia8)

install(
  TARGETS ActsExamplesGeneratorsPythia8
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
