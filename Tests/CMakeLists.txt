# helpers shared between all tests
add_subdirectory(CommonHelpers)

# add a unittest executable w/ default dependencies and register it
macro(add_unittest _name)
  # for now we use the same name also for the target
  set(_target ${_name})
  add_executable(${_target} ${ARGN})
  target_include_directories(
    ${_target}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(
    ${_target}
    PRIVATE ActsCore ActsTestsCommonHelpers Boost::unit_test_framework)
  # register as unittest executable
  add_test(NAME ${_name} COMMAND ${_target})
endmacro()

add_subdirectory_if(Core ACTS_BUILD_TESTS)
add_subdirectory_if(Examples ACTS_BUILD_EXAMPLES)
add_subdirectory_if(IntegrationTests ACTS_BUILD_INTEGRATION_TESTS)
add_subdirectory_if(Plugins ACTS_BUILD_TESTS)
add_subdirectory_if(Benchmarks ACTS_BUILD_BENCHMARKS)
