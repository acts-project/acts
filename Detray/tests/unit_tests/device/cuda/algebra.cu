// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Test include(s).
#include "detray/test/device/cuda/algebra_test_suite.cuh"

// GoogleTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <string>

namespace detray::test::cuda {

using test_types = testing::Types<detray::test::algebra>;

REGISTER_TYPED_TEST_SUITE_P(cuda_vector_test, vector_2d_ops, vector_3d_ops);
REGISTER_TYPED_TEST_SUITE_P(cuda_matrix_test, matrix64_ops, matrix22_ops);
REGISTER_TYPED_TEST_SUITE_P(cuda_transform_test, transform3D);

// The 'test_types' are defined in 'algebra/test/framework/types.hpp'
INSTANTIATE_TYPED_TEST_SUITE_P(detray_algebra, cuda_vector_test, test_types,
                               test_specialisation_name);
INSTANTIATE_TYPED_TEST_SUITE_P(detray_algebra, cuda_matrix_test, test_types,
                               test_specialisation_name);
INSTANTIATE_TYPED_TEST_SUITE_P(detray_algebra, cuda_transform_test, test_types,
                               test_specialisation_name);

}  // namespace detray::test::cuda
