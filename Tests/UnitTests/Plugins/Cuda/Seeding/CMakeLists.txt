add_executable(SeedfinderCudaTest SeedfinderCudaTest.cpp)
target_link_libraries(SeedfinderCudaTest PRIVATE ${unittest_extra_libraries} Boost::boost)

