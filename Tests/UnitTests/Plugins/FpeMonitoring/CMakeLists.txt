set(unittest_extra_libraries ActsPluginFpeMonitoring)

if(NOT APPLE)
  # This test doesn't make sense on macOS at this time
  add_unittest(FpeMonitor FpeMonitorTests.cpp)
  if(("${ACTS_BUILD_ONLY_MATCHING_UNITTESTS}" STREQUAL "") OR ("ActsUnitTestFpeMonitor" MATCHES "${ACTS_BUILD_ONLY_MATCHING_UNITTESTS}"))
    set_tests_properties(FpeMonitor PROPERTIES ENVIRONMENT BOOST_TEST_CATCH_SYSTEM_ERRORS=no)
  endif()
endif()
