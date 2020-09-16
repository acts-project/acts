
add_executable(ActsUnitTestSeedfinderCuda2
  main.cpp
  CommandLineArguments.hpp CommandLineArguments.cpp
  ReadSeedFile.hpp ReadSeedFile.cpp
  TestSpacePoint.hpp TestSpacePoint.cpp
  TestDeviceCuts.hpp TestDeviceCuts.cu
  TestHostCuts.hpp TestHostCuts.cpp)
target_link_libraries(ActsUnitTestSeedfinderCuda2
  PRIVATE Boost::program_options ActsCore ActsPluginCuda2)
set_target_properties(ActsUnitTestSeedfinderCuda2
  PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
