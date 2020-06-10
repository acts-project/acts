add_executable(SeedFinderCudaTest SeedFinderCudaTest.cpp)
target_link_libraries(SeedFinderCudaTest PRIVATE ${unittest_extra_libraries} Boost::boost)

