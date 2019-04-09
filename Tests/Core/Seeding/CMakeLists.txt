add_executable (SeedfinderTest SeedfinderTest.cpp)
target_link_libraries (SeedfinderTest PRIVATE ActsCore)
add_test (NAME SeedfinderTests COMMAND SeedfinderTest)
acts_add_test_to_cdash_project (PROJECT ACore TEST SeedfinderTests TARGETS SeedfinderTest)
