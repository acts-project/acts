set(unittest_extra_libraries ActsPluginMille)

# Need the `pede` program for the solver test
if(NOT ACTS_USE_SYSTEM_MP2)
    find_program(
        PEDE_EXECUTABLE
        NAMES pede
        HINTS ${millepedeII_BINARY_DIR}/bin
        DOC "Millepede II solver binary"
    )
endif()

add_unittest(ActsKalmanToMilleConsistency ActsKalmanToMilleConsistency.cpp)
add_unittest(ActsMilleBasicCalls ActsMilleBasicCalls.cpp)
add_unittest(MillePedeSolverBuiltIn MillePedeSolverBuiltIn.cpp)

target_compile_definitions(
    ActsUnitTestMillePedeSolverBuiltIn
    PRIVATE
        ACTS_PEDE_EXECUTABLE="${PEDE_EXECUTABLE}";
        ACTS_LOG_FAILURE_THRESHOLD=FATAL
)
add_unittest(ChildProcessTests ChildProcessTests.cpp)
