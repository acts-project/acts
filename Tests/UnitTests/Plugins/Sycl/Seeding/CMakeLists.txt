add_executable(ActsUnitTestSeedfinderSycl
  SeedfinderSyclTest.cpp
  CommandLineArguments.h
  CommandLineArguments.cpp
  ATLASCuts.hpp
  SpacePoint.hpp)
target_link_libraries(ActsUnitTestSeedfinderSycl
  PRIVATE Boost::program_options ActsCore ActsPluginSycl)
  