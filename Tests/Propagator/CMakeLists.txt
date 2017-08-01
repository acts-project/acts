add_executable (PropagationTests PropagationTests.cpp)
target_link_libraries (PropagationTests PRIVATE ACTSCore)
add_test(NAME PropagationTests COMMAND PropagationTests)

add_executable (RungeKuttaEngineTests RungeKuttaEngineTests.cpp)
target_link_libraries (RungeKuttaEngineTests PRIVATE ACTSCore)
add_test(NAME RungeKuttaEngineTests COMMAND RungeKuttaEngineTests)

