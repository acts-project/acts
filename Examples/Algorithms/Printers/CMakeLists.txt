add_library(
  ActsExamplesPrinters SHARED
  ACTFW/Printers/PrintHits.cpp
  ACTFW/Printers/PrintParticles.cpp)
target_include_directories(
  ActsExamplesPrinters
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(
  ActsExamplesPrinters
  PUBLIC ActsCore ActsFatras ActsDigitizationPlugin ActsExamplesFramework)

install(
  TARGETS ActsExamplesPrinters
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
