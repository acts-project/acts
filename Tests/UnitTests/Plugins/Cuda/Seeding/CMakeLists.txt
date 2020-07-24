
# The unit test needs to explicitly use the CUDA Runtime API.
find_library(CUDART_LIBRARY cudart ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})

add_executable(SeedfinderCudaTest SeedfinderCudaTest.cpp)
target_link_libraries(SeedfinderCudaTest
  PRIVATE ActsCore ActsCudaPlugin Boost::boost ${CUDART_LIBRARY})
target_include_directories(SeedfinderCudaTest
  PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
