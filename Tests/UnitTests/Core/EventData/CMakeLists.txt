add_unittest(BoundTrackParameters BoundTrackParametersTests.cpp)
add_unittest(Charge ChargeTests.cpp)
add_unittest(CurvilinearTrackParameters CurvilinearTrackParametersTests.cpp)
add_unittest(FreeTrackParameters FreeTrackParametersTests.cpp)
add_unittest(MeasurementHelpers MeasurementHelpersTests.cpp)
add_unittest(Measurement MeasurementTests.cpp)
add_unittest(MultiComponentBoundTrackParameters MultiComponentBoundTrackParametersTests.cpp)
add_unittest(MultiTrajectory MultiTrajectoryTests.cpp)
add_unittest(TransformBoundToFree TransformBoundToFreeTests.cpp)
add_unittest(TransformFreeToBound TransformFreeToBoundTests.cpp)
add_unittest(CorrectedTransformFreeToBound CorrectedTransformFreeToBoundTests.cpp)
add_unittest(Track TrackTests.cpp)
add_unittest(SourceLink SourceLinkTests.cpp)

function(add_non_compile_test name testname src)
  set(target "${name}${testname}")
  add_executable(${target} ${src})
  set_target_properties(${target} PROPERTIES
    EXCLUDE_FROM_ALL ON 
    EXCLUDE_FROM_DEFAULT_BUILD ON)
  target_compile_definitions(${target} PRIVATE _TEST="${testname}")
  target_link_libraries(${target} PUBLIC ActsCore)

  add_test(NAME ActsNonCompileTest${target}
    COMMAND ${CMAKE_COMMAND} --build . --target ${target} --config $<CONFIGURATION>
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  set_tests_properties(ActsNonCompileTest${target} PROPERTIES WILL_FAIL ON)
endfunction()

set(src "${CMAKE_CURRENT_SOURCE_DIR}/MultiTrajectoryComplianceTests.cpp")

file(READ ${src} content)

string(REGEX MATCHALL "ACTS_DOES_NOT_COMPILE_BEGIN\\(([A-Za-z0-9]+)\\)" matches ${content})
foreach(match ${matches})
  string(REGEX REPLACE "ACTS_DOES_NOT_COMPILE_BEGIN\\(([A-Za-z0-9]+)\\)" "\\1" match ${match})
  if(${match} STREQUAL "x")
    continue()
  endif()
  message(STATUS ${match})
  add_non_compile_test(MultiTrajectory ${match} ${src})
endforeach()


# add_non_compile_test(MultiTrajectory Closure ${src})
# set_tests_properties(ActsNonCompileTestClosure PROPERTIES WILL_FAIL ON)

set(name "MultiTrajectory")
set(target "${name}Closure")
add_executable(${target} ${src})
set_target_properties(${target} PROPERTIES
  EXCLUDE_FROM_ALL ON 
  EXCLUDE_FROM_DEFAULT_BUILD ON)
target_compile_definitions(${target} PRIVATE _TEST="XXXXNOTEST__")
target_link_libraries(${target} PUBLIC ActsCore)

add_test(NAME ActsNonCompileTest${target}
  COMMAND ${CMAKE_COMMAND} --build . --target ${target} --config $<CONFIGURATION>
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(ActsNonCompileTest${target} PROPERTIES WILL_FAIL OFF)
