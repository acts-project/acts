# shared code 
add_subdirectory(Common)

# Digitization for Generic detector
add_executable(
  ActsExampleDigitizationConfigGeneric
  GenericDigitizationConfigExample.cpp)
target_link_libraries(
  ActsExampleDigitizationConfigGeneric
  PRIVATE
    ActsExamplesCommon
    ActsExamplesDigitizationCommon
)

install(
  TARGETS
  ActsExampleDigitizationConfigGeneric 
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(
  ActsExampleDigitizationGeneric
  GenericDigitizationExample.cpp)
target_link_libraries(
  ActsExampleDigitizationGeneric
  PRIVATE
    ActsExamplesCommon
    ActsExamplesDigitizationCommon
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
    ActsExamplesCommon
    ActsExamplesDigitizationCommon
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
    ActsExamplesCommon
    ActsExamplesDigitizationCommon
)

install(
  TARGETS
    ActsExampleDigitizationTGeo 
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_subdirectory_if(DD4hep ACTS_BUILD_EXAMPLES_DD4HEP)
