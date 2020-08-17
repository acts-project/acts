add_library(
  ActsExamplesGenerators SHARED
  ActsExamples/Generators/EventGenerator.cpp
  ActsExamples/Generators/FlattenEvent.cpp
  ActsExamples/Generators/ParametricProcessGenerator.cpp
  ActsExamples/Generators/ParticleSelector.cpp
  ActsExamples/Options/ParticleGunOptions.cpp)
target_include_directories(
  ActsExamplesGenerators
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(
  ActsExamplesGenerators
  PUBLIC ActsCore ActsExamplesFramework Boost::program_options)

install(
  TARGETS ActsExamplesGenerators
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
