add_executable(SeedfinderTest SeedfinderTest.cpp)
target_link_libraries(SeedfinderTest PRIVATE ActsCore Boost::boost)

if (ACTS_BUILD_CUDA_FEATURES)  
add_executable(SeedfinderCUDAValidate SeedfinderCUDAValidate.cpp)
target_link_libraries(SeedfinderCUDAValidate PRIVATE ${CUDART_LIBRARY} ActsCore Boost::boost)
endif()
