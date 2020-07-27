
add_executable(SeedfinderCudaTest
  SeedfinderCudaTest.hpp SeedfinderCudaTest.cpp SeedfinderCudaTest.cu
  ATLASCuts.hpp SpacePoint.hpp)
target_link_libraries(SeedfinderCudaTest
  PRIVATE ActsCore ActsCudaPlugin)
