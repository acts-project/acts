
# prepare integration tests
# add custom target that runs the invokation shell script
option (BUILD_INTEGRATION_TESTS "build integration tests" OFF)

if (BUILD_INTEGRATION_TESTS)

  add_custom_target(integration_tests
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
    )

  # macro to collect integration test binaries from subdirectories
  macro(add_integration_test)
    file (RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
    add_custom_target("Integration${ARGV0}" 
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      COMMAND "${_relPath}/${ARGV0}")
    add_dependencies(integration_tests "Integration${ARGV0}")
  endmacro()
  
  
  add_executable (PropagationTests PropagationTests.cpp)
  target_link_libraries (PropagationTests PRIVATE ACTSCore)
  add_integration_test(PropagationTests)

else (BUILD_INTEGRATION_TESTS)
  message (STATUS "${Blue}disable build of Integration Tests${ColorReset}")
endif (BUILD_INTEGRATION_TESTS)
