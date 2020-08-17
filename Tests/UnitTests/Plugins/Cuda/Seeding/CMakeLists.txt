add_executable(ActsUnitTestSeedfinderCuda SeedfinderCudaTest.cpp)
target_link_libraries(ActsUnitTestSeedfinderCuda PRIVATE ${unittest_extra_libraries} Boost::boost)

