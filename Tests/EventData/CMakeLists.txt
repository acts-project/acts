# allow to turn of unit tests manually
option (BUILD_TESTS "build unit tests" ON)

if (BUILD_TESTS)
   # we need boost
   find_package (Boost REQUIRED)

   if (Boost_FOUND)
      add_executable (ParameterSetTest ParameterSetTests.cpp)
      target_include_directories (ParameterSetTest PRIVATE ${Boost_INCLUDE_DIRS})
      target_link_libraries (ParameterSetTest PRIVATE ACTSCore)
      add_test(NAME ParSetUnitTest COMMAND ParameterSetTest)

      add_executable (MeasurementTest MeasurementTests.cpp)
      target_include_directories (MeasurementTest PRIVATE ${Boost_INCLUDE_DIRS})
      target_link_libraries (MeasurementTest PRIVATE ACTSCore)
      add_test(NAME MeasurementUnitTest COMMAND MeasurementTest)
   endif (Boost_FOUND)
else (BUILD_TESTS)
  message (STATUS "${Blue}disable build unit tests{ColorReset}")
endif (BUILD_TESTS)
