add_library(
  ActsCustomLogger SHARED
  CustomDefaultLogger.cpp)
target_include_directories(
  ActsCustomLogger
  PRIVATE ${PROJECT_SOURCE_DIR}/Core/include)

add_executable(
  ActsTabulateEnergyLoss
  TabulateEnergyLoss.cpp)
target_link_libraries(
  ActsTabulateEnergyLoss
  PRIVATE ActsCore ActsTestsCommonHelpers)

install(
  TARGETS ActsCustomLogger ActsTabulateEnergyLoss
  EXPORT ActsExamplesTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  COMPONENT Examples)
