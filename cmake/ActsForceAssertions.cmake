# Force assertions regardless of build type when ACTS_FORCE_ASSERTIONS is set.
#
# This is shared by ACTS and the in-tree Detray and Traccc projects, so that
# their standalone builds honour the same option. It only affects targets
# created in the including directory and below, so it has to be included
# before any targets are set up. The include guard makes nested inclusions a
# no-op: a project added via add_subdirectory() inherits the settings from the
# directory that included this first.
include_guard(GLOBAL)

if(ACTS_FORCE_ASSERTIONS)
    message(
        STATUS
        "Injecting headers to force assertions. This can have side-effects, USE WITH CAUTION!"
    )
    include_directories(BEFORE SYSTEM ${CMAKE_CURRENT_LIST_DIR}/assert_include)

    # Add assertions to standard libraries. These go into the flags rather than
    # the directory compile definitions, which ROOT_GENERATE_DICTIONARY copies
    # into a custom command where a $<COMPILE_LANGUAGE> restriction is invalid.
    string(APPEND CMAKE_CXX_FLAGS " -D_GLIBCXX_ASSERTIONS -D_LIBCPP_DEBUG")
endif()

# Register death tests that check the forced assertions actually fire. Each
# probe is built with NDEBUG defined and must fail with an assertion message.
#
#   acts_add_assertion_tests([CUDA])
#
# CUDA adds a probe for device-side assert(). Every probe is added only once per
# build, so ACTS, Detray and Traccc can all call this when built together.
function(acts_add_assertion_tests)
    if(NOT ACTS_FORCE_ASSERTIONS)
        return()
    endif()

    cmake_parse_arguments(PARSE_ARGV 0 _args "CUDA" "" "")

    set(_probes AssertProbe.cpp StdlibAssertProbe.cpp)
    if(_args_CUDA)
        list(APPEND _probes CudaAssertProbe.cu)
    endif()

    foreach(_source IN LISTS _probes)
        get_filename_component(_name ${_source} NAME_WE)
        get_property(_added GLOBAL PROPERTY _ACTS_ASSERTION_TEST_${_name})
        if(_added)
            continue()
        endif()
        set_property(GLOBAL PROPERTY _ACTS_ASSERTION_TEST_${_name} TRUE)

        set(_target ActsForceAssertions${_name})
        add_executable(
            ${_target}
            ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/assert_check/${_source}
        )
        target_compile_definitions(${_target} PRIVATE NDEBUG)
        add_test(
            NAME ForceAssertions.${_name}
            COMMAND
                ${CMAKE_COMMAND} -DPROBE=$<TARGET_FILE:${_target}> -P
                ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/assert_check/RunAssertProbe.cmake
        )
        set_tests_properties(
            ForceAssertions.${_name}
            PROPERTIES SKIP_REGULAR_EXPRESSION "ACTS_ASSERTION_PROBE_SKIPPED"
        )
    endforeach()
endfunction()
