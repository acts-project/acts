
add_executable(SeedFinderCuda2Test
  main.cpp
  CommandLineArguments.hpp CommandLineArguments.cpp
  ReadSeedFile.hpp ReadSeedFile.cpp
  TestSpacePoint.hpp TestSpacePoint.cpp
  TestDeviceCuts.hpp TestDeviceCuts.cu
  TestHostCuts.hpp TestHostCuts.cpp)
target_link_libraries(SeedFinderCuda2Test
  PRIVATE Boost::program_options ActsCore ActsCuda2Plugin)
set_target_properties(SeedFinderCuda2Test
  PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
