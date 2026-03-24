set(unittest_extra_libraries ActsPluginFpeMonitoring)

add_unittest(FpeMonitor FpeMonitorTests.cpp)
set_tests_properties(
    FpeMonitor
    PROPERTIES ENVIRONMENT BOOST_TEST_CATCH_SYSTEM_ERRORS=no
)
