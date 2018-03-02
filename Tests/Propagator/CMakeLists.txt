add_executable (ExtendableTests ExtendableTests.cpp)
target_link_libraries (ExtendableTests PRIVATE ACTSCore)
add_test(NAME ExtendableUnitTests COMMAND ExtendableTests)

add_executable (AbortListTests AbortListTests.cpp)
target_link_libraries (AbortListTests PRIVATE ACTSCore)
add_test(NAME AbortListUnitTest COMMAND AbortListTests)

add_executable (ActionListTests ActionListTests.cpp)
target_link_libraries (ActionListTests PRIVATE ACTSCore)
add_test(NAME ActionListUnitTests COMMAND ActionListTests)

add_executable (JacobianTests JacobianTests.cpp)
target_link_libraries (JacobianTests PRIVATE ACTSCore)
add_test(NAME JacobianUnitTests COMMAND JacobianTests)

add_executable (PropagatorTests PropagatorTests.cpp)
target_link_libraries (PropagatorTests PRIVATE ACTSCore)
add_test(NAME PropagatorUnitTests COMMAND PropagatorTests)

add_executable (RungeKuttaEngineTests RungeKuttaEngineTests.cpp)
target_link_libraries (RungeKuttaEngineTests PRIVATE ACTSCore)
add_test(NAME RungeKuttaEngineUnitTest COMMAND RungeKuttaEngineTests)

# Benchmark tests

# define executables
add_executable (AtlasStepperBenchmark AtlasStepperBenchmark.cpp)
target_include_directories (AtlasStepperBenchmark PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (AtlasStepperBenchmark PRIVATE ACTSCore)
target_link_libraries (AtlasStepperBenchmark PRIVATE ${Boost_LIBRARIES})

add_executable (EigenStepperBenchmark EigenStepperBenchmark.cpp)
target_include_directories (EigenStepperBenchmark PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (EigenStepperBenchmark PRIVATE ACTSCore)
target_link_libraries (EigenStepperBenchmark PRIVATE ${Boost_LIBRARIES})

add_executable (AtlasPropagatorBenchmark AtlasPropagatorBenchmark.cpp)
target_include_directories (AtlasPropagatorBenchmark PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (AtlasPropagatorBenchmark PRIVATE ACTSCore)
target_link_libraries (AtlasPropagatorBenchmark PRIVATE ${Boost_LIBRARIES})

install (TARGETS AtlasStepperBenchmark EigenStepperBenchmark AtlasPropagatorBenchmark
	     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/Profiling)
