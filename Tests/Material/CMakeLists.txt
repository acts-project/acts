add_executable (MaterialTests MaterialTests.cpp)
target_include_directories (MaterialTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (MaterialTests PRIVATE ACTSCore)
add_test (NAME MaterialUnitTest COMMAND MaterialTests)
acts_add_test_to_cdash_project (PROJECT ACore TEST MaterialUnitTest TARGETS MaterialTests)

