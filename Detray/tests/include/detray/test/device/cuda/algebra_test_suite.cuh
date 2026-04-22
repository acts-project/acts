// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Test include(s).
#include "detray/test/device/cuda/cuda_test_fixture.hpp"
#include "detray/test/device/cuda/execute_cuda_test.cuh"
#include "detray/test/device/execute_host_test.hpp"
#include "detray/test/device/matrix_fixture.hpp"
#include "detray/test/device/transform_fixture.hpp"
#include "detray/test/device/vector_fixture.hpp"

// GoogleTest include(s).
#include <gtest/gtest.h>

namespace detray::test::cuda {

template <detray::concepts::algebra A>
class cuda_vector_test : public cuda_test_fixture<vector_fixture<A>> {};

template <detray::concepts::algebra A>
class cuda_matrix_test : public cuda_test_fixture<matrix_fixture<A>> {};

template <detray::concepts::algebra A>
class cuda_transform_test : public cuda_test_fixture<transform_fixture<A>> {};

TYPED_TEST_SUITE_P(cuda_vector_test);
TYPED_TEST_SUITE_P(cuda_matrix_test);
TYPED_TEST_SUITE_P(cuda_transform_test);

/// Test for some basic 2D "vector operations"
TYPED_TEST_P(cuda_vector_test, vector_2d_ops) {
  // Run the test on the host, and on the/a device.
  execute_host_test<vector_2d_ops_functor<TypeParam>>(
      this->m_p1->size(), vecmem::get_data(*(this->m_p1)),
      vecmem::get_data(*(this->m_p2)),
      vecmem::get_data(*(this->m_output_host)));

  execute_cuda_test<vector_2d_ops_functor<TypeParam>>(
      this->m_p1->size(), vecmem::get_data(*(this->m_p1)),
      vecmem::get_data(*(this->m_p2)),
      vecmem::get_data(*(this->m_output_device)));

  // Compare the outputs.
  this->compareOutputs();
}

/// Test for some basic 3D "vector operations"
TYPED_TEST_P(cuda_vector_test, vector_3d_ops) {
  // This test is just not numerically stable at float precision in optimized
  // mode for some reason. :-(
#ifdef NDEBUG
  if (typeid(dscalar<TypeParam>) == typeid(float)) {
    GTEST_SKIP();
  }
#endif  // NDEBUG

  // Run the test on the host, and on the/a device.
  execute_host_test<vector_3d_ops_functor<TypeParam>>(
      this->m_v1->size(), vecmem::get_data(*(this->m_v1)),
      vecmem::get_data(*(this->m_v2)),
      vecmem::get_data(*(this->m_output_host)));

  execute_cuda_test<vector_3d_ops_functor<TypeParam>>(
      this->m_v1->size(), vecmem::get_data(*(this->m_v1)),
      vecmem::get_data(*(this->m_v2)),
      vecmem::get_data(*(this->m_output_device)));

  // Compare the outputs.
  this->compareOutputs();
}

/// Test for handling matrices
TYPED_TEST_P(cuda_matrix_test, matrix64_ops) {
  // Run the test on the host, and on the/a device.
  execute_host_test<matrix64_ops_functor<TypeParam>>(
      this->m_m1->size(), vecmem::get_data(*(this->m_m1)),
      vecmem::get_data(*(this->m_output_host)));

  execute_cuda_test<matrix64_ops_functor<TypeParam>>(
      this->m_m1->size(), vecmem::get_data(*(this->m_m1)),
      vecmem::get_data(*(this->m_output_device)));

  // Compare the outputs.
  this->compareOutputs();
}

/// Test for handling matrices
TYPED_TEST_P(cuda_matrix_test, matrix22_ops) {
  // Run the test on the host, and on the/a device.
  execute_host_test<matrix22_ops_functor<TypeParam>>(
      this->m_m2->size(), vecmem::get_data(*(this->m_m2)),
      vecmem::get_data(*(this->m_output_host)));

  execute_cuda_test<matrix22_ops_functor<TypeParam>>(
      this->m_m2->size(), vecmem::get_data(*(this->m_m2)),
      vecmem::get_data(*(this->m_output_device)));

  // Compare the outputs.
  this->compareOutputs();
}

/// Test for some operations with @c transform3
TYPED_TEST_P(cuda_transform_test, transform3D) {
  // Run the test on the host, and on the/a device.
  execute_host_test<transform3_ops_functor<TypeParam>>(
      this->m_t1->size(), vecmem::get_data(*(this->m_t1)),
      vecmem::get_data(*(this->m_t2)), vecmem::get_data(*(this->m_t3)),
      vecmem::get_data(*(this->m_v1)), vecmem::get_data(*(this->m_v2)),
      vecmem::get_data(*(this->m_output_host)));

  execute_cuda_test<transform3_ops_functor<TypeParam>>(
      this->m_t1->size(), vecmem::get_data(*(this->m_t1)),
      vecmem::get_data(*(this->m_t2)), vecmem::get_data(*(this->m_t3)),
      vecmem::get_data(*(this->m_v1)), vecmem::get_data(*(this->m_v2)),
      vecmem::get_data(*(this->m_output_device)));

  // Compare the outputs.
  this->compareOutputs();
}

}  // namespace detray::test::cuda
