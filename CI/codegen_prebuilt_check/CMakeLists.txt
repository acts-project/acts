# Checks that a directory of pre-generated code is complete and that
# acts_code_generation actually picks it up.
#
# This deliberately is not ACTS: it declares no languages and builds no code, so
# it needs neither a compiler nor any of the C++ dependencies a real ACTS
# configure requires. What it does exercise is the real ActsCodegen.cmake, the
# real manifest and the real lookup, with ACTS_CODEGEN_REQUIRE_PREBUILT turning
# any fall back to running a generator into a configure error.
#
#   CI/pregenerate_codegen.py --output <prebuilt>
#   cmake -S CI/codegen_prebuilt_check -B <build> \
#         -DACTS_SOURCE_DIR=<acts> -DACTS_CODEGEN_PREBUILT_DIR=<prebuilt>
#   cmake --build <build>
#
# Leaving ACTS_CODEGEN_PREBUILT_DIR unset instead checks that the directory is
# found on its own next to the top-level CMakeLists.txt, which is how a build
# from the release source archive finds it.

cmake_minimum_required(VERSION 3.19)
project(ActsCodegenPrebuiltCheck NONE)

# This project lives at <acts>/CI/codegen_prebuilt_check, so the source tree is
# known without being told. ACTS_SOURCE_DIR is only there to point the check at
# a different tree, such as an unpacked source archive; pass it as an absolute
# path, since a relative one resolves against this directory.
if(NOT DEFINED ACTS_SOURCE_DIR)
    set(ACTS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../..")
endif()
get_filename_component(ACTS_SOURCE_DIR "${ACTS_SOURCE_DIR}" ABSOLUTE)

if(NOT EXISTS "${ACTS_SOURCE_DIR}/cmake/ActsCodegen.cmake")
    message(
        FATAL_ERROR
        "ACTS_SOURCE_DIR=${ACTS_SOURCE_DIR} does not look like an ACTS source tree"
    )
endif()

set(ACTS_CODEGEN_REQUIRE_PREBUILT ON CACHE BOOL "" FORCE)
include(${ACTS_SOURCE_DIR}/cmake/ActsCodegen.cmake)

acts_codegen_manifest(_manifest)
string(JSON _count LENGTH "${_manifest}" units)
if(_count EQUAL 0)
    message(FATAL_ERROR "The codegen manifest is empty")
endif()
math(EXPR _last "${_count} - 1")

foreach(_i RANGE ${_last})
    string(JSON _key MEMBER "${_manifest}" units ${_i})

    # Stand-in for the real consumer: acts_code_generation only asks its target
    # for the TYPE and then attaches an include directory and a dependency, so
    # an INTERFACE library exercises the same wiring ActsCore does.
    string(MAKE_C_IDENTIFIER "${_key}" _target)
    add_library(${_target} INTERFACE)

    acts_code_generation(ADD_TO_TARGET ${_target} KEY ${_key})
endforeach()

message(STATUS "Checked ${_count} codegen units")
