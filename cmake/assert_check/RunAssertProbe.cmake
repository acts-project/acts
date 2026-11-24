# Run an assertion probe and check that it died on a failed assertion. CTest
# reports a test killed by a signal as failed regardless of its output, so the
# probe cannot be registered directly.
#
#   cmake -DPROBE=<executable> -P RunAssertProbe.cmake

execute_process(
    COMMAND ${PROBE}
    RESULT_VARIABLE _result
    OUTPUT_VARIABLE _output
    ERROR_VARIABLE _output
)
message("${_output}")

if("${_result}" STREQUAL "77")
    message("ACTS_ASSERTION_PROBE_SKIPPED")
elseif("${_result}" STREQUAL "0")
    message(FATAL_ERROR "${PROBE} exited normally, assertions are not enabled")
elseif(NOT _output MATCHES "Assertion.*failed")
    message(
        FATAL_ERROR
        "${PROBE} failed (${_result}) without an assertion message"
    )
endif()
