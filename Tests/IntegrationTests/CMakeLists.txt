# integration tests are different and should not be registered with the
# preconfigured `tests` target. manually add a separate one.
add_custom_target(integrationtests)
# add an integrationtest exectuable w/ default dependencies and register it.
#
# the common libraries which are linked to every integrationtest can be
# extended by setting the `unittest_extra_libraries` variables before
# calling the macro.
macro(add_integrationtest _name)
  # automatically prefix the target name
  set(_target "ActsIntegrationTest${_name}")
  add_executable(${_target} ${ARGN})
  # define required BOOST_TEST_... macros here to ensure consistent names
  target_compile_definitions(
    ${_target}
    PRIVATE "-DBOOST_TEST_DYN_LINK" "-DBOOST_TEST_MODULE=${_target}")
  target_include_directories(
    ${_target}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(
    ${_target}
    PRIVATE
      ActsCore
      ActsTestsCommonHelpers
      Boost::unit_test_framework
      ${integrationtest_extra_libraries})
  # register as integration test executable
  set(_test "integrationtests-${_name}")
  add_custom_target(${_test} COMMAND ${_target})
  add_dependencies(integrationtests ${_test})
endmacro()

add_integrationtest(BVHNavigation BVHNavigationTest.cpp)
add_integrationtest(InterpolatedSolenoidBField InterpolatedSolenoidBFieldTest.cpp)
add_integrationtest(Propagation PropagationTests.cpp)

add_subdirectory_if(Fatras ACTS_BUILD_FATRAS)
add_subdirectory_if(Legacy ACTS_BUILD_LEGACY)
