# Digitization for DD4hep detector
add_executable(
  ActsExampleDigitizationConfigDD4hep
  DD4hepDigitizationConfigExample.cpp)
target_link_libraries(
  ActsExampleDigitizationConfigDD4hep
  PRIVATE
    ActsExamplesDigitizationCommon
    ActsExamplesDD4hepDetectorWithOptions
)

install(
  TARGETS
  ActsExampleDigitizationConfigDD4hep 
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(
  ActsExampleDigitizationDD4hep
  DD4hepDigitizationExample.cpp)
target_link_libraries(
    ActsExampleDigitizationDD4hep
  PRIVATE
    ActsExamplesDigitizationCommon
    ActsExamplesDD4hepDetectorWithOptions)

install(
  TARGETS ActsExampleDigitizationDD4hep
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
