find_package (Boost COMPONENTS unit_test_framework REQUIRED)
find_package (Eigen3 REQUIRED)

include_directories (${Boost_INCLUDE_DIRS}
		     ${EIGEN_INCLUDE_DIR}
		     ${CMAKE_SOURCE_DIR}/Core/Algebra
		     ${CMAKE_SOURCE_DIR}/Core/CoreIdentifier
		     ${CMAKE_SOURCE_DIR}/Core/CoreUtils
		     ${CMAKE_SOURCE_DIR}/EventData/ParameterSet
		     ${CMAKE_SOURCE_DIR}/EventData/EventDataUtils
		     ${CMAKE_SOURCE_DIR}/EventData/ParametersBase
		     ${CMAKE_SOURCE_DIR}/EventData/TrackParameters
		     ${CMAKE_SOURCE_DIR}/EventData/Measurement
		     ${CMAKE_SOURCE_DIR}/Geometry/Surfaces
		     ${CMAKE_SOURCE_DIR}/Geometry/DetectorElementBase
		     ${CMAKE_SOURCE_DIR}/Geometry/GeometryUtils
		    )

link_directories(${PROJECT_BINARY_DIR}/lib
				/afs/cern.ch/exp/fcc/sw/0.5/LCG_80/tbb/42_20140122/x86_64-slc6-gcc49-opt/lib
				/afs/cern.ch/exp/fcc/sw/0.5/LCG_80/Boost/1.59.0_python2.7/x86_64-slc6-gcc49-opt/lib
				/afs/cern.ch/exp/fcc/sw/0.5/GAUDI/GAUDI_v26r4/InstallArea/x86_64-slc6-gcc49-opt/lib
				)
add_definitions (-DBOOST_TEST_DYN_LINK)
add_executable (ParameterSetTest ParameterSetTests.cxx)
add_executable (MeasurementTest MeasurementTests.cxx)
target_link_libraries (ParameterSetTest ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
target_link_libraries (MeasurementTest ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} libSurfaces.so)
add_test(NAME ParSetUnitTest COMMAND ParameterSetTest)
add_test(NAME MeasurementUnitTest COMMAND MeasurementTest)
