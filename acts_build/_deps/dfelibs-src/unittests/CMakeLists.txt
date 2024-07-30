# minimal version that provides test v3
find_package(Boost 1.59 REQUIRED COMPONENTS unit_test_framework)
# optional, for namedtuple_python test
find_package(PythonInterp 2.7)
# optional, for io_root test
find_package(ROOT 6.10)

function(add_unittest _name)
  set(_target "${PROJECT_NAME}_unittest_${_name}")
  add_executable(${_target} "test_${_name}.cpp")
  target_compile_definitions(${_target} PRIVATE DFE_UNITTESTS_DIRECTORY="${CMAKE_CURRENT_SOURCE_DIR}")
  target_compile_definitions(${_target} PRIVATE BOOST_TEST_MODULE="${_name}")
  target_compile_definitions(${_target} PRIVATE BOOST_TEST_DYN_LINK=1)
  # as system include to suppress warnings
  target_include_directories(${_target} SYSTEM PRIVATE ${Boost_INCLUDE_DIR})
  target_link_libraries(${_target} PRIVATE dfelibs)
  target_link_libraries(${_target} PRIVATE ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
  add_test(NAME ${_name} COMMAND ${_target})
endfunction()

add_unittest(dispatcher)
add_unittest(flatmap)
add_unittest(flatset)
add_unittest(histogram)
add_unittest(io_dsv)
add_unittest(io_numpy)
if(ROOT_FOUND)
  add_unittest(io_root)
  # ROOT might require C++17 but does not advertise it
  set(_target "${PROJECT_NAME}_unittest_io_root")
  target_compile_features(${_target} PRIVATE cxx_std_17)
  target_link_libraries(${_target} PRIVATE ROOT::Tree)
endif()
if(PYTHONINTERP_FOUND)
  add_test(
    NAME io_read_from_python
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_io_read_from_python.py)
  # python test use the output written by the namedtuple test
  set_tests_properties(io_dsv PROPERTIES FIXTURES_SETUP io_read_from_python)
  set_tests_properties(io_numpy PROPERTIES FIXTURES_SETUP io_read_from_python)
  set_tests_properties(io_read_from_python PROPERTIES FIXTURES_REQUIRED io_read_from_python)
endif()
add_unittest(namedtuple)

add_unittest(poly)
add_unittest(smallflatset)
add_unittest(smallvector)
