add_executable(SeedFinderSyclTest SeedFinderSyclTest.cpp)
target_link_libraries(SeedFinderSyclTest PRIVATE ${unittest_extra_libraries} Boost::boost)