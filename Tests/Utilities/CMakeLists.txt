add_executable (unit_conversion_tests unit_conversion.cpp)
target_link_libraries (unit_conversion_tests PRIVATE ACTSCore)
add_test(NAME unit_conversion_tests COMMAND unit_conversion_tests)
add_dependencies(ACoreProject unit_conversion_tests)

add_executable (GeometryID_tests GeometryID_tests.cpp)
target_link_libraries (GeometryID_tests PRIVATE ACTSCore)
add_test(NAME GeometryID_tests COMMAND GeometryID_tests)
add_dependencies(ACoreProject GeometryID_tests)

add_executable (BinningData_tests BinningData_tests.cpp)
target_link_libraries (BinningData_tests PRIVATE ACTSCore)
add_test(NAME BinningData_tests COMMAND BinningData_tests)
add_dependencies(ACoreProject BinningData_tests)

add_executable (BinUtility_tests BinUtility_tests.cpp)
target_link_libraries (BinUtility_tests PRIVATE ACTSCore)
add_test(NAME BinUtility_tests COMMAND BinUtility_tests)
add_dependencies(ACoreProject BinUtility_tests)

add_executable (Logger_tests Logger_tests.cpp)
target_link_libraries (Logger_tests PRIVATE ACTSCore)
add_test(NAME Logger_tests COMMAND Logger_tests)
add_dependencies(ACoreProject Logger_tests)
