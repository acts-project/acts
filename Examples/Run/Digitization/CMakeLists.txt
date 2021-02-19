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
