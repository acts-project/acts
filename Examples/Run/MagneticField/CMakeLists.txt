add_executable(
  ActsExampleMagneticField
  BFieldExample.cpp)
target_link_libraries(
  ActsExampleMagneticField
  PRIVATE
    ActsCore
    ActsExamplesFramework ActsExamplesCommon
    ActsExamplesMagneticField ActsExamplesIoRoot
    Boost::program_options)

add_executable(
  ActsExampleMagneticFieldAcess
  BFieldAccessExample.cpp)
target_link_libraries(
  ActsExampleMagneticFieldAcess
  PRIVATE
    ActsCore
    ActsExamplesFramework ActsExamplesCommon
    ActsExamplesMagneticField ActsExamplesIoRoot Boost::program_options)

install(
  TARGETS ActsExampleMagneticField ActsExampleMagneticFieldAcess
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
