add_library(
  ActsExamplesGenerators SHARED
  ACTFW/Generators/EventGenerator.cpp
  ACTFW/Generators/FlattenEvent.cpp
  ACTFW/Generators/ParametricProcessGenerator.cpp
  ACTFW/Generators/ParticleSelector.cpp
  ACTFW/Options/ParticleGunOptions.cpp)
target_include_directories(
  ActsExamplesGenerators
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(
  ActsExamplesGenerators
  PUBLIC ActsCore ActsExamplesFramework Boost::program_options)

install(
  TARGETS ActsExamplesGenerators
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
