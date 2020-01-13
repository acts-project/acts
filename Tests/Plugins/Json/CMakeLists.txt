add_executable (JsonConversionTest JsonConversionTest.cpp)
target_include_directories (JsonConversionTest PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (JsonConversionTest PRIVATE ActsCore ActsJsonPlugin ActsTestsCommonHelpers ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test (NAME JsonConversionTest COMMAND JsonConversionTest)
