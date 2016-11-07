add_executable (ParameterSetTest ParameterSetTests.cpp)
target_include_directories (ParameterSetTest PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (ParameterSetTest PRIVATE ACore)
add_test(NAME ParSetUnitTest COMMAND ParameterSetTest)

add_executable (MeasurementTest MeasurementTests.cpp)
target_include_directories (MeasurementTest PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (MeasurementTest PRIVATE ACore)
add_test(NAME MeasurementUnitTest COMMAND MeasurementTest)
