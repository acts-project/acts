add_executable(
  ActsTabulateEnergyLoss
  TabulateEnergyLoss.cpp)
target_link_libraries(
  ActsTabulateEnergyLoss
  PRIVATE ActsCore ActsFatras)

install(
  TARGETS ActsTabulateEnergyLoss
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
