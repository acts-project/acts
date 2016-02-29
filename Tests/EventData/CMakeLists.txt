find_package (Boost COMPONENTS unit_test_framework REQUIRED)
find_package (Eigen3 REQUIRED)

include_directories (${Boost_INCLUDE_DIRS}
		     ${EIGEN_INCLUDE_DIR}
		     ${CMAKE_SOURCE_DIR}/Core/Algebra
		     ${CMAKE_SOURCE_DIR}/Core/CoreUtils
		     ${CMAKE_SOURCE_DIR}/EventData/ParameterSet
		    )

add_definitions (-DBOOST_TEST_DYN_LINK)
add_executable (ParameterSetTest ParameterSetTests.cxx)
target_link_libraries (ParameterSetTest ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME ParSetUnitTest COMMAND ParameterSetTest)
