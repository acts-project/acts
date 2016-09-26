add_executable (SeedingTests SeedingTests.cpp)
target_include_directories (SeedingTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (SeedingTests PRIVATE ACTSCore)
add_test(NAME SeedingUnitTests COMMAND SeedingTests)

add_executable (SeedingToolsTests SeedingToolsTests.cpp)
target_include_directories (SeedingToolsTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (SeedingToolsTests PRIVATE ACTSCore)
add_test(NAME SeedingToolsUnitTests COMMAND SeedingToolsTests)
