add_executable (ExtendableTests ExtendableTests.cpp)
target_link_libraries (ExtendableTests PRIVATE ActsCore ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME ExtendableUnitTests COMMAND ExtendableTests)

add_executable (UnitConversionTests UnitConversionTests.cpp)
target_link_libraries (UnitConversionTests PRIVATE ActsCore ActsTestsCommonHelpers ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test (NAME UnitConversionUnitTest COMMAND UnitConversionTests)
acts_add_test_to_cdash_project (PROJECT ACore TEST UnitConversionUnitTest TARGETS UnitConversionTests)

add_executable (GeometryIDTests GeometryIDTests.cpp)
target_link_libraries (GeometryIDTests PRIVATE ActsCore ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test (NAME GeometryIDUnitTest COMMAND GeometryIDTests)
acts_add_test_to_cdash_project (PROJECT ACore TEST GeometryIDUnitTest TARGETS GeometryIDTests)

add_executable (IntersectionTests IntersectionTests.cpp)
target_link_libraries (IntersectionTests PRIVATE ActsCore ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test (NAME IntersectionUnitTests COMMAND IntersectionTests)
acts_add_test_to_cdash_project (PROJECT ACore TEST IntersectionUnitTests TARGETS IntersectionTests)

add_executable (BinningDataTests BinningDataTests.cpp)
target_link_libraries (BinningDataTests PRIVATE ActsCore ActsTestsCommonHelpers ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test (NAME BinningDataUnitTest COMMAND BinningDataTests)
acts_add_test_to_cdash_project (PROJECT ACore TEST BinningDataUnitTest TARGETS BinningDataTests)

add_executable (BinUtilityTests BinUtilityTests.cpp)
target_link_libraries (BinUtilityTests PRIVATE ActsCore ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test (NAME BinUtilityUnitTest COMMAND BinUtilityTests)
acts_add_test_to_cdash_project (PROJECT ACore TEST BinUtilityUnitTest TARGETS BinUtilityTests)

add_executable (LoggerTests LoggerTests.cpp)
target_link_libraries (LoggerTests PRIVATE ActsCore ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test (NAME LoggerUnitTest COMMAND LoggerTests)
acts_add_test_to_cdash_project (PROJECT ACore TEST LoggerUnitTest TARGETS LoggerTests)

add_executable (MPLTests MPLTests.cpp)
target_link_libraries (MPLTests PRIVATE ActsCore ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test (NAME MPLUnitTest COMMAND MPLTests)
acts_add_test_to_cdash_project (PROJECT ACore TEST MPLUnitTest TARGETS MPLTests)

add_executable (RealQuadraticEquationTests RealQuadraticEquationTests.cpp)
target_link_libraries (RealQuadraticEquationTests PRIVATE ActsCore ActsTestsCommonHelpers ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test (NAME RealQuadraticEquationUnitTest COMMAND RealQuadraticEquationTests)

add_executable (AxesTests AxesTests.cpp)
target_include_directories (AxesTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (AxesTests PRIVATE ActsCore ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test (NAME AxesUnitTest COMMAND AxesTests)
acts_add_test_to_cdash_project (PROJECT ACore TEST AxesUnitTest TARGETS AxesTests)

add_executable (GridTests GridTests.cpp)
target_include_directories (GridTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (GridTests PRIVATE ActsCore ActsTestsCommonHelpers ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test (NAME GridUnitTest COMMAND GridTests)
acts_add_test_to_cdash_project (PROJECT ACore TEST GridUnitTest TARGETS GridTests)

add_executable (InterpolationTests InterpolationTests.cpp)
target_include_directories (InterpolationTests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (InterpolationTests PRIVATE ActsCore ActsTestsCommonHelpers ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test (NAME InterpolationUnitTest COMMAND InterpolationTests)
acts_add_test_to_cdash_project (PROJECT ACore TEST InterpolationUnitTest TARGETS InterpolationTests)

add_executable (BFieldMapUtilsTests BFieldMapUtilsTests.cpp)
target_link_libraries (BFieldMapUtilsTests PRIVATE ActsCore ActsTestsCommonHelpers ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test (NAME BFieldMapUtilsTest COMMAND BFieldMapUtilsTests)
acts_add_test_to_cdash_project (PROJECT ACore TEST BFieldMapUtilsTest TARGETS BFieldMapUtilsTests)

add_executable (ResultTests ResultTests.cpp)
target_link_libraries (ResultTests PRIVATE ActsCore ActsTestsCommonHelpers ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test (NAME ResultTests COMMAND ResultTests)

add_executable (TypeTraitsTest TypeTraitsTest.cpp)
target_link_libraries (TypeTraitsTest PRIVATE ActsCore ActsTestsCommonHelpers ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test (NAME TypeTraitsTest COMMAND TypeTraitsTest)

add_executable (VisualizationTests VisualizationTests.cpp)
target_link_libraries (VisualizationTests PRIVATE ActsCore ActsTestsCommonHelpers ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test (NAME VisualizationTests COMMAND VisualizationTests)

add_executable (MaterialMapUtilsTests MaterialMapUtilsTests.cpp)
target_link_libraries (MaterialMapUtilsTests PRIVATE ActsCore ActsTestsCommonHelpers)
add_test (NAME MaterialMapUtilsTest COMMAND MaterialMapUtilsTests)
acts_add_test_to_cdash_project (PROJECT ACore TEST MaterialMapUtilsTest TARGETS MaterialMapUtilsTests)