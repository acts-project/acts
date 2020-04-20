add_executable(SeedfinderCudaValidate SeedfinderCudaValidate.cpp)
target_link_libraries(SeedfinderCudaValidate PRIVATE ${unittest_extra_libraries} ActsCore Boost::boost ${CUDART_LIBRARY})

