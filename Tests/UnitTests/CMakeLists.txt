# add a unittest executable w/ default dependencies and register it
macro(add_unittest _name)
  # automatically prefix the target name
  set(_target "ActsUnitTest${_name}")
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

add_subdirectory(Core)
add_subdirectory(Plugins)
