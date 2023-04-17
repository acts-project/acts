# add a unittest executable w/ default dependencies and register it.

# the common libraries which are linked to every unittest can be
# extended by setting the `unittest_extra_libraries` variables before
# calling the macro.
macro(add_unittest _name)
  # automatically prefix the target name
  set(_target "ActsUnitTest${_name}")
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
      ${unittest_extra_libraries})
  # register as unittest executable
  add_test(NAME ${_name} COMMAND ${_target})
endmacro()

function(add_non_compile_test name src)

  if(NOT ACTS_BUILD_NONCOMPILE_TESTS)
    return()
  endif()

  cmake_path(ABSOLUTE_PATH src)
  get_filename_component(_filename ${src} NAME)
  set(_processed_source "${CMAKE_CURRENT_BINARY_DIR}/${_filename}")
  
  add_custom_command(
    OUTPUT ${_processed_source}
    DEPENDS ${src}
    COMMAND "${CMAKE_COMMAND}"
            -DINPUT_FILE=${src}
            -DOUTPUT_FILE=${_processed_source}
            -P ${CMAKE_SOURCE_DIR}/cmake/ActsGenerateNonCompileTest.cmake
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
  add_custom_target(
      "${name}_generated_source"
      DEPENDS ${_processed_source}
  )

  set(test "ActsNonCompileTest${name}Closure")
  set(target "${test}_Executable")
  add_executable(${target} ${_processed_source})
  target_link_libraries(${target} PUBLIC ActsCore ActsTestsCommonHelpers)
  add_dependencies(${target} "${name}_generated_source")

  set_target_properties(${target} PROPERTIES
    EXCLUDE_FROM_ALL ON 
    EXCLUDE_FROM_DEFAULT_BUILD ON)

  add_test(NAME ${test}
    COMMAND ${CMAKE_COMMAND} --build . --target ${target} --config $<CONFIGURATION>
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  file(READ ${src} content)
  string(REGEX MATCHALL "ACTS_DOES_NOT_COMPILE_BEGIN\\(([A-Za-z0-9]+)\\)" matches ${content})
  foreach(match ${matches})
    string(REGEX REPLACE "ACTS_DOES_NOT_COMPILE_BEGIN\\(([A-Za-z0-9]+)\\)" "\\1" match ${match})

    set(test "ActsNonCompileTest${name}${match}")
    set(target "${test}_Executable")
    add_executable(${target} ${_processed_source})
    target_link_libraries(${target} PUBLIC ActsCore ActsTestsCommonHelpers)
    target_compile_definitions(${target} PRIVATE "-D${match}")
    add_dependencies(${target} "${name}_generated_source")

    set_target_properties(${target} PROPERTIES
      EXCLUDE_FROM_ALL ON 
      EXCLUDE_FROM_DEFAULT_BUILD ON)

    add_test(NAME ${test}
      COMMAND ${CMAKE_COMMAND} --build . --target ${target} --config $<CONFIGURATION>
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
    set_tests_properties(${test} PROPERTIES WILL_FAIL ON)
  endforeach()

endfunction()

add_subdirectory(Core)
add_subdirectory_if(Examples ACTS_BUILD_EXAMPLES)
add_subdirectory_if(Benchmarks ACTS_BUILD_BENCHMARKS)
add_subdirectory_if(Fatras ACTS_BUILD_FATRAS)
add_subdirectory(Plugins)
add_subdirectory_if(Alignment ACTS_BUILD_ALIGNMENT)
