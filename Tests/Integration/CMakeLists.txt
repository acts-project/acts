add_custom_target(integration_tests WORKING_DIRECTORY ${CMAKE_PROJECT_DIR})

# macro to collect integration test binaries from subdirectories
macro(add_integration_test)
  file (RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
  add_custom_target("Integration${ARGV0}"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMAND "${_relPath}/${ARGV0}")
  add_dependencies(integration_tests "Integration${ARGV0}")
endmacro()

add_executable (PropagationTests PropagationTests.cpp)
target_link_libraries (PropagationTests PRIVATE ActsCore ActsLegacy ActsTestsCommonHelpers)
add_integration_test(PropagationTests)

add_executable (SeedingTest ATLSeedingIntegrationTest.cpp)
target_link_libraries (SeedingTest PRIVATE ActsLegacy)
add_integration_test(SeedingTest)

add_executable (InterpolatedSolenoidBFieldTest InterpolatedSolenoidBFieldTest.cpp)
target_link_libraries (InterpolatedSolenoidBFieldTest PRIVATE ActsCore ActsTestsCommonHelpers)
add_integration_test(InterpolatedSolenoidBFieldTest)
