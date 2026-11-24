acts_apply_warning_flags(CXX ${ACTS_FAIL_ON_WARNINGS})

# helpers shared between all tests
add_subdirectory(CommonHelpers)

add_subdirectory_if(Benchmarks ACTS_BUILD_BENCHMARKS)
add_subdirectory_if(IntegrationTests ACTS_BUILD_INTEGRATIONTESTS)
add_subdirectory_if(UnitTests ACTS_BUILD_UNITTESTS)
