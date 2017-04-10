add_executable (unit_conversion_tests unit_conversion.cpp)
target_link_libraries (unit_conversion_tests PRIVATE ACTSCore)
add_test (NAME unit_conversion_tests COMMAND unit_conversion_tests)
acts_add_test_to_cdash_project (PROJECT ACore TEST unit_conversion_tests TARGETS unit_conversion_tests)

add_executable (GeometryID_tests GeometryID_tests.cpp)
target_link_libraries (GeometryID_tests PRIVATE ACTSCore)
add_test (NAME GeometryID_tests COMMAND GeometryID_tests)
acts_add_test_to_cdash_project (PROJECT ACore TEST GeometryID_tests TARGETS GeometryID_tests)

add_executable (BinningData_tests BinningData_tests.cpp)
target_link_libraries (BinningData_tests PRIVATE ACTSCore)
add_test (NAME BinningData_tests COMMAND BinningData_tests)
acts_add_test_to_cdash_project (PROJECT ACore TEST BinningData_tests TARGETS BinningData_tests)

add_executable (BinUtility_tests BinUtility_tests.cpp)
target_link_libraries (BinUtility_tests PRIVATE ACTSCore)
add_test (NAME BinUtility_tests COMMAND BinUtility_tests)
acts_add_test_to_cdash_project (PROJECT ACore TEST BinUtility_tests TARGETS BinUtility_tests)

add_executable (Logger_tests Logger_tests.cpp)
target_link_libraries (Logger_tests PRIVATE ACTSCore)
add_test (NAME Logger_tests COMMAND Logger_tests)
acts_add_test_to_cdash_project (PROJECT ACore TEST Logger_tests TARGETS Logger_tests)

add_executable (MPL_tests MPL_checks.cpp)
target_link_libraries (MPL_tests PRIVATE ACTSCore)
add_test (NAME MPL_tests COMMAND MPL_tests)
acts_add_test_to_cdash_project (PROJECT ACore TEST MPL_tests TARGETS MPL_tests)

add_executable (RealQuadraticEquationTests RealQuadraticEquationTests.cpp)
target_link_libraries (RealQuadraticEquationTests PRIVATE ACTSCore)
add_test (NAME RealQuadraticEquationTests COMMAND RealQuadraticEquationTests)
