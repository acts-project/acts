add_library(
  ActsExamplesPrinters SHARED
  ActsExamples/Printers/PrintHits.cpp
  ActsExamples/Printers/PrintParticles.cpp)
target_include_directories(
  ActsExamplesPrinters
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(
  ActsExamplesPrinters
  PUBLIC ActsCore ActsFatras ActsPluginDigitization ActsExamplesFramework)

install(
  TARGETS ActsExamplesPrinters
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
