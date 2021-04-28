add_executable(
  ActsExampleDigitizationDD4hep
  DD4hepDigitizationExample.cpp)
target_link_libraries(
    ActsExampleDigitizationDD4hep
  PRIVATE ActsExamplesDigitizationCommon ActsExamplesDetectorDD4hep)

install(
  TARGETS ActsExampleDigitizationDD4hep
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
