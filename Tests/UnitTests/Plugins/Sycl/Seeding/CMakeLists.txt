add_executable(ActsUnitTestSeedfinderSycl SeedfinderSyclTest.cpp CommandLineArguments.cpp)
FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(ActsUnitTestSeedfinderSycl PRIVATE ${unittest_extra_libraries} ${Boost_LIBRARIES})