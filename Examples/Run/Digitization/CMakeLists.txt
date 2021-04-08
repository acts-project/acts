# shared code 
add_subdirectory(Common)

# Digitization for Generic detector
add_executable(
  ActsExampleDigitizationGeneric
  GenericDigitizationExample.cpp)
target_link_libraries(
  ActsExampleDigitizationGeneric
  PRIVATE
    ActsExamplesDigitizationCommon
    ActsExamplesDetectorGeneric
)

install(
  TARGETS
    ActsExampleDigitizationGeneric
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# Digitization for TGeo detector
add_executable(
  ActsExampleDigitizationConfigTGeo
  TGeoDigitizationConfigExample.cpp)
target_link_libraries(
  ActsExampleDigitizationConfigTGeo
  PRIVATE
    ActsExamplesDigitizationCommon
    ActsExamplesDetectorTGeo
)

install(
  TARGETS
    ActsExampleDigitizationConfigTGeo 
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(
  ActsExampleDigitizationTGeo
  TGeoDigitizationExample.cpp)
target_link_libraries(
  ActsExampleDigitizationTGeo
  PRIVATE
    ActsExamplesDigitizationCommon
    ActsExamplesDetectorTGeo
)

install(
  TARGETS
    ActsExampleDigitizationTGeo 
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})