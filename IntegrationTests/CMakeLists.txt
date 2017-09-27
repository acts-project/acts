
# prepare integration tests
# add custom target that runs the invokation shell script
option (BUILD_INTEGRATION_TESTS "build integration tests" OFF)

if (BUILD_INTEGRATION_TESTS)

  find_package (ROOT ${ACTS_ROOT_VERSION} ${ACTS_FIND_QUIET_TEXT} COMPONENTS RIO Tree MathCore)
  if (ROOT_FOUND)


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

  else (ROOT_FOUND)
    message (STATUS "${Red}ROOT was not found -> skip build of Integration Tests${ColorReset}")
  endif (ROOT_FOUND)

else (BUILD_INTEGRATION_TESTS)
  message (STATUS "${Blue}disable build of Integration Tests${ColorReset}")
endif (BUILD_INTEGRATION_TESTS)
