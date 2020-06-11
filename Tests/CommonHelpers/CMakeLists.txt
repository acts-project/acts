# get the path of the test data directory
file(TO_NATIVE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../Data" acts_test_data_dir)

add_library(
  ActsTestsCommonHelpers SHARED
  Acts/Tests/CommonHelpers/DataDirectory.cpp)
target_compile_definitions(
  ActsTestsCommonHelpers
  PRIVATE "ACTS_TEST_DATA_DIR=\"${acts_test_data_dir}\"")
target_compile_features(
  ActsTestsCommonHelpers
  PUBLIC cxx_std_17)
target_include_directories(
  ActsTestsCommonHelpers
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(
  ActsTestsCommonHelpers
  INTERFACE ActsCore
  PRIVATE std::filesystem)
