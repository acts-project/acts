# SPDX-License-Identifier: MIT
# Copyright 2020 Moritz Kiehn
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

# Provide `add_format_targets` to add `format` and `check-format` targets.
#
# To use the module, include it in your CMakeLists.txt and call
# `add_format_targets` with all source files/ source globs that should be
# formatted
#
#     include(ClangFormatTargets)
#     add_format_targets(src/file1.cpp src/*.hpp)
#
# The `add_format_targets` is always available, but can be a no-op if
# clang-format could not be found. The `format` target will be created, if
# clang-format is available. The `check-format` target also requires git.

find_program(CLANG_FORMAT_EXECUTABLE clang-format)
find_package(Git)

# ensure the provided function is always available
if(NOT CLANG_FORMAT_EXECUTABLE)
  function(add_format_targets)
  endfunction()
  return()
endif()

message(STATUS "Found clang-format: ${CLANG_FORMAT}")

function(add_format_targets)
  file(GLOB_RECURSE sources ${ARGN})
  list(SORT sources)

  message(DEBUG "Added format sources:")
  foreach(source ${sources})
    message(DEBUG "  ${source}")
  endforeach()

  add_custom_target(
    format
    COMMAND ${CLANG_FORMAT_EXECUTABLE} --style=file -i ${sources}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Formatting files files with clang-format"
    SOURCES ${sources})

  # check target requires git for comparison with the checked-in files
  if(GIT_EXECUTABLE)
    add_custom_target(
      check-format
      COMMAND ${GIT_EXECUTABLE} diff --stat --exit-code
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Check format compliance")
    add_dependencies(check-format format)
  endif()
endfunction()
