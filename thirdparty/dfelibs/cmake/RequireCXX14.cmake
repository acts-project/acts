# SPDX-License-Identifier: MIT
# Copyright 2015-2020 Moritz Kiehn
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

# Require support for C++14 globally.
#
# If a higher standard has already been set via CMAKE_CXX_FLAGS it is retained,
# existing flags for lower standards are removed.
#
# WARNING
# Starting with CMake 3.8, this module is deprecated. A much cleaner way to
# require C++14 (or any other language standard) is with per-target compile
# features. Use e.g.
#
#     target_compile_features(<target> INTERFACE cxx_std_14)
#

include(CheckCXXCompilerFlag)

check_cxx_compiler_flag("-std=c++14" _has_cxx14_support)
if(NOT _has_cxx14_support)
  message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++14 support.")
endif()

# check for existing compiler flags
string(REGEX MATCHALL "-std=[a-z0-9+]+" _defined_stds ${CMAKE_CXX_FLAGS})
set(_compatible_stds)
foreach(_std ${_defined_stds})
  # extract standard version number
  string(REGEX REPLACE "-std=(c|gnu)\\+\\+([a-z0-9]+)" "\\2" _version ${_std})
  if(_version MATCHES "98|03|0x|11|1y")
    message(STATUS "Remove C++ standard flag ${_std}")
  else()
    message(STATUS "Keep C++ standard flag ${_std}")
    set(_compatible_stds "${_compatible_stds} ${_std}")
  endif()
endforeach()

if(_compatible_stds)
  # (some) existing standard flags are compatible w/ c++14
  # remove all standard flags
  string(REGEX REPLACE "-std=[a-z0-9+]+" "" _cxx_flags ${CMAKE_CXX_FLAGS})
  # readd compatible ones
  set(CMAKE_CXX_FLAGS "${_cxx_flags} ${_compatible_stds}")
else()
  # existing standard flags are not enough
  if(CMAKE_VERSION VERSION_LESS "3.1")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
    message(STATUS "Added -std=c++14 compiler flag.")
  else()
    set(CMAKE_CXX_STANDARD 14)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
    message(STATUS "Require CXX_STANDARD=14.")
  endif()
endif()
