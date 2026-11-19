# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

message(STATUS "Building the std::array SoA algebra plugin")

# Set up the library.
detray_add_library( detray_algebra_array_soa algebra_array_soa
   "include/algebra/array_soa.hpp"
   "include/algebra/impl/*.hpp"
   "include/algebra/utils/*.hpp"
)

target_link_libraries(detray_algebra_array_soa INTERFACE detray::core)

# The lane-wise math functions only vectorize on the host, if they do not have
# to set errno. Detray does not read errno, so the flag does not change the
# results. Device code never sets errno, so the device compilation needs no flag.
set(_detray_gnu_like_cxx
    "$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:IntelLLVM>>"
)
target_compile_options(
    detray_algebra_array_soa
    INTERFACE
        "$<$<AND:$<COMPILE_LANGUAGE:CXX>,${_detray_gnu_like_cxx}>:-fno-math-errno>"
        "$<$<AND:$<COMPILE_LANGUAGE:CUDA>,${_detray_gnu_like_cxx}>:SHELL:-Xcompiler -fno-math-errno>"
        "$<$<COMPILE_LANGUAGE:HIP>:-fno-math-errno>"
)
unset(_detray_gnu_like_cxx)

target_compile_definitions(
    detray_algebra_array_soa
    INTERFACE
        DETRAY_CUSTOM_SCALARTYPE=${DETRAY_CUSTOM_SCALARTYPE}
        DETRAY_ALGEBRA_ARRAY_SOA
)

# Compile all headers one-by-one.
detray_compile_headers(
    algebra_array_soa
    detray::algebra_array_soa
    "include/algebra/*.hpp"
    "include/algebra/**/*.hpp"
)
