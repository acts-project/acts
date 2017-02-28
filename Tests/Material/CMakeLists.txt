add_executable (Material_tests Material_tests.cpp)
target_include_directories (Material_tests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (Material_tests PRIVATE ACTSCore)
add_test (NAME Material_tests COMMAND Material_tests)
acts_add_test_to_cdash_project (PROJECT ACore TEST Material_tests TARGETS Material_tests)

