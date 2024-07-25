# optional, for namedtuple example
find_package(ROOT 6.10)

function(add_example _name)
  set(_target "${PROJECT_NAME}_example_${_name}")
  add_executable(${_target} "example_${_name}.cpp")
  target_link_libraries(${_target} PRIVATE dfelibs)
endfunction()

add_example(dispatcher)
add_example(namedtuple)
if (ROOT_FOUND)
  set(_target "${PROJECT_NAME}_example_namedtuple")
  target_compile_definitions(${_target} PRIVATE DFE_USE_IO_ROOT)
  # ROOT might require C++17 but does not advertise it
  target_compile_features(${_target} PRIVATE cxx_std_17)
  target_link_libraries(${_target} PRIVATE ROOT::Tree)
endif()
