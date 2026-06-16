# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

message(STATUS "Building detray SYCL integration tests")

# Enable SYCL as a language.
include(vecmem-check-language) # Needed for the SYCL language.
enable_language(SYCL)

# Set the SYCL build flags.
include(detray-compiler-options-sycl)

# make unit tests for multiple algebras
# Currently eigen, vc and smatrix is not supported
set(algebras "array")

foreach(algebra ${algebras})
    # Unit tests for the selected algebra.
    detray_add_integration_test(sycl_${algebra}
       "propagator_sycl_kernel.hpp"
       "propagator.sycl"
       "propagator_kernel.sycl"
       LINK_LIBRARIES GTest::gtest_main vecmem::sycl detray::test_device
       detray::core detray::algebra_${algebra} detray::test_utils
    )
    target_compile_definitions(
        detray_integration_test_sycl_${algebra}
        PRIVATE ${algebra}=${algebra}
    )
endforeach()
