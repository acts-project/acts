add_library(
  ActsExampleCustomLogger SHARED
  CustomDefaultLogger.cpp)
target_compile_features(
  ActsExampleCustomLogger
  PUBLIC cxx_std_17)

target_include_directories(
  ActsExampleCustomLogger
  PRIVATE ${PROJECT_SOURCE_DIR}/Core/include)

add_executable(
  ActsTabulateEnergyLoss
  TabulateEnergyLoss.cpp)
target_link_libraries(
  ActsTabulateEnergyLoss
  PRIVATE ActsCore ActsFatras)

install(
  TARGETS ActsExampleCustomLogger ActsTabulateEnergyLoss
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
