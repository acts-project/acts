add_executable(SeedFinderSyclTest SeedFinderSyclTest.cpp CommandLineArguments.cpp)
FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(SeedFinderSyclTest PRIVATE ${unittest_extra_libraries} ${Boost_LIBRARIES})