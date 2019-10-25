add_library(ActsTestsCommonHelpers INTERFACE)
target_include_directories(
  ActsTestsCommonHelpers
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(
  ActsTestsCommonHelpers
  INTERFACE ActsCore)
