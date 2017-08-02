add_executable (PropagationTests PropagationTests.cpp)
target_link_libraries (PropagationTests PRIVATE ACTSCore)
add_test(NAME PropagationUnitTest COMMAND PropagationTests)

add_executable (RungeKuttaEngineTests RungeKuttaEngineTests.cpp)
target_link_libraries (RungeKuttaEngineTests PRIVATE ACTSCore)
add_test(NAME RungeKuttaEngineUnitTest COMMAND RungeKuttaEngineTests)

# Benchmark tests

# we need boost
find_package (Boost ${ACTS_BOOST_VERSION} REQUIRED COMPONENTS program_options)

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


