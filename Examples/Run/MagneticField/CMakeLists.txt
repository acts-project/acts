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
  ActsExampleMagneticFieldAccess
  BFieldAccessExample.cpp)
target_link_libraries(
  ActsExampleMagneticFieldAccess
  PRIVATE
    ActsCore
    ActsExamplesFramework ActsExamplesCommon
    ActsExamplesMagneticField ActsExamplesIoRoot)

install(
  TARGETS ActsExampleMagneticField ActsExampleMagneticFieldAccess
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
