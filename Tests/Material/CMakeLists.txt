add_executable (MaterialTests MaterialTests.cpp)
target_include_directories (MaterialTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (MaterialTests PRIVATE ActsCore)

add_executable (MaterialPropertiesTests MaterialPropertiesTests.cpp)
target_include_directories (MaterialPropertiesTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (MaterialPropertiesTests PRIVATE ActsCore)

add_test (NAME MaterialUnitTest COMMAND MaterialTests MaterialPropertiesTests)
acts_add_test_to_cdash_project (PROJECT ACore TEST MaterialUnitTest TARGETS MaterialTests MaterialPropertiesTests)
