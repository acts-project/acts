set(unittest_extra_libraries ActsPluginMille)

# Need the `pede` program for the solver test.
# This is not found when running the "ACTS-internal" MP2 build
# as part of the ACTS CI. Hence this test is only run
# when we have a complete Millepede-II installation in place.
if(ACTS_USE_SYSTEM_MP2)
    add_unittest(MillePedeSolverBuiltIn MillePedeSolverBuiltIn.cpp)
endif()

add_unittest(ActsKalmanToMilleConsistency ActsKalmanToMilleConsistency.cpp)
add_unittest(ActsMilleBasicCalls ActsMilleBasicCalls.cpp)
add_unittest(ChildProcessTests ChildProcessTests.cpp)
