add_library(
  ActsExamplesDigitizationCommon SHARED
  DigitizationConfigExample.cpp
  DigitizationExample.cpp
  DigitizationInput.cpp)
target_include_directories(
  ActsExamplesDigitizationCommon
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(
  ActsExamplesDigitizationCommon
  PUBLIC
    ActsExamplesFramework
    ActsExamplesCommon
    ActsExamplesDigitization
    ActsExamplesIoCsv
    ActsExamplesIoJson
    ActsExamplesIoRoot)

install(
  TARGETS ActsExamplesDigitizationCommon
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
