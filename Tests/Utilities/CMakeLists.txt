add_executable (unit_conversion_tests unit_conversion.cpp)
add_executable (GeometryID_tests GeometryID_tests.cpp)
add_executable (BinningData_tests BinningData_tests.cpp)
add_executable (BinUtility_tests BinUtility_tests.cpp)
target_link_libraries (unit_conversion_tests PRIVATE ACTSCore)
target_link_libraries (GeometryID_tests PRIVATE ACTSCore)
target_link_libraries (BinningData_tests PRIVATE ACTSCore)
target_link_libraries (BinUtility_tests PRIVATE ACTSCore)
add_test(NAME unit_conversion_tests COMMAND unit_conversion_tests)
add_test(NAME GeometryID_tests COMMAND GeometryID_tests)
add_test(NAME BinUtility_tests COMMAND BinUtility_tests)


