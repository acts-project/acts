add_executable(SeedfinderTest SeedfinderTest.cpp)
add_definitions(-Wno-implicit-fallthrough)
target_link_libraries(SeedfinderTest PRIVATE ActsCore Boost::boost)
