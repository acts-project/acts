add_library(
  ActsExamplesFatrasCommon SHARED
  src/FatrasCommon.cpp src/FatrasOptions.cpp)
target_include_directories(
  ActsExamplesFatrasCommon
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesFatrasCommon
  PUBLIC
    ActsCore
    ActsExamplesFramework ActsExamplesCommon
    ActsExamplesGenerators
    ActsExamplesMagneticField
    ActsExamplesFatras
    ActsExamplesIoCsv ActsExamplesIoRoot)

install(
  TARGETS ActsExamplesFatrasCommon
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
