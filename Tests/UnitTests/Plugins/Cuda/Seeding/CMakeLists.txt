add_executable(SeedfinderCudaValidate SeedfinderCudaValidate.cpp)
target_link_libraries(SeedfinderCudaValidate PRIVATE ${unittest_extra_libraries} Boost::boost)

