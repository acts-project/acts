# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

message(STATUS "Building detray SYCL unit tests")

# Enable SYCL as a language.
include(vecmem-check-language) # Needed for the SYCL language.
enable_language(SYCL)

# Set the default SYCL compiler flags.
include(detray-compiler-options-sycl)

# Define unit tests for multiple algebra plugins
# Currently, only 'array' is supported in SYCL code
set(algebras "array")

foreach(algebra ${algebras})
    # Unit tests for the selected algebra.
    detray_add_unit_test(sycl_${algebra}
       "algebra.sycl"
       LINK_LIBRARIES GTest::gtest_main vecmem::sycl detray::test_common detray::test_sycl
    )

    target_compile_definitions(
        detray_unit_test_sycl_${algebra}
        PRIVATE ${algebra}=${algebra}
    )
endforeach()
