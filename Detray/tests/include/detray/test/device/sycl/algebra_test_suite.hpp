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
#include "detray/test/device/execute_host_test.hpp"
#include "detray/test/device/matrix_fixture.hpp"
#include "detray/test/device/sycl/execute_sycl_test.hpp"
#include "detray/test/device/sycl/sycl_test_fixture.hpp"
#include "detray/test/device/transform_fixture.hpp"
#include "detray/test/device/vector_fixture.hpp"

// GoogleTest include(s).
#include <gtest/gtest.h>

// SYCL include(s).
#include <sycl/sycl.hpp>

namespace detray::test::sycl {

/// Test case class, to be specialised for the different plugins
template <detray::concepts::algebra A>
class sycl_vector_test : public sycl_test_fixture<vector_fixture<A>> {};

/// Test case class, to be specialised for the different plugins
template <detray::concepts::algebra A>
class sycl_matrix_test : public sycl_test_fixture<matrix_fixture<A>> {};

/// Test case class, to be specialised for the different plugins
template <detray::concepts::algebra A>
class sycl_transform_test : public sycl_test_fixture<transform_fixture<A>> {};

TYPED_TEST_SUITE_P(sycl_vector_test);
TYPED_TEST_SUITE_P(sycl_matrix_test);
TYPED_TEST_SUITE_P(sycl_transform_test);

/// Test for some basic 2D "vector operations"
TYPED_TEST_P(sycl_vector_test, vector_2d_ops) {
  // Don't run the test at double precision, if the SYCL device doesn't
  // support it.
  if ((typeid(dscalar<TypeParam>) == typeid(double)) &&
      (this->m_queue.get_device().has(::sycl::aspect::fp64) == false)) {
    GTEST_SKIP();
  }

  // Run the test on the host, and on the/a device.
  execute_host_test<vector_2d_ops_functor<TypeParam>>(
      this->m_p1->size(), vecmem::get_data(*(this->m_p1)),
      vecmem::get_data(*(this->m_p2)),
      vecmem::get_data(*(this->m_output_host)));

  execute_sycl_test<vector_2d_ops_functor<TypeParam>>(
      this->m_queue, this->m_p1->size(), vecmem::get_data(*(this->m_p1)),
      vecmem::get_data(*(this->m_p2)),
      vecmem::get_data(*(this->m_output_device)));

  // Compare the outputs.
  this->compareOutputs();
}

/// Test for some basic 3D "vector operations"
TYPED_TEST_P(sycl_vector_test, vector_3d_ops) {
  // Don't run the test at double precision, if the SYCL device doesn't
  // support it.
  if ((typeid(dscalar<TypeParam>) == typeid(double)) &&
      (this->m_queue.get_device().has(::sycl::aspect::fp64) == false)) {
    GTEST_SKIP();
  }

  // This test is just not numerically stable at float precision in optimized
  // mode on some backends. :-( (Cough... HIP... cough...)
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

  execute_sycl_test<vector_3d_ops_functor<TypeParam>>(
      this->m_queue, this->m_v1->size(), vecmem::get_data(*(this->m_v1)),
      vecmem::get_data(*(this->m_v2)),
      vecmem::get_data(*(this->m_output_device)));

  // Compare the outputs.
  this->compareOutputs();
}
/// Test for handling matrices
TYPED_TEST_P(sycl_matrix_test, matrix64_ops) {
  // Don't run the test at double precision, if the SYCL device doesn't
  // support it.
  if ((typeid(dscalar<TypeParam>) == typeid(double)) &&
      (this->m_queue.get_device().has(::sycl::aspect::fp64) == false)) {
    GTEST_SKIP();
  }

  // Run the test on the host, and on the/a device.
  execute_host_test<matrix64_ops_functor<TypeParam>>(
      this->m_m1->size(), vecmem::get_data(*(this->m_m1)),
      vecmem::get_data(*(this->m_output_host)));

  execute_sycl_test<matrix64_ops_functor<TypeParam>>(
      this->m_queue, this->m_m1->size(), vecmem::get_data(*(this->m_m1)),
      vecmem::get_data(*(this->m_output_device)));

  // Compare the outputs.
  this->compareOutputs();
}

/// Test for handling matrices
TYPED_TEST_P(sycl_matrix_test, matrix22_ops) {
  // Don't run the test at double precision, if the SYCL device doesn't
  // support it.
  if ((typeid(dscalar<TypeParam>) == typeid(double)) &&
      (this->m_queue.get_device().has(::sycl::aspect::fp64) == false)) {
    GTEST_SKIP();
  }

  // Run the test on the host, and on the/a device.
  execute_host_test<matrix22_ops_functor<TypeParam>>(
      this->m_m2->size(), vecmem::get_data(*(this->m_m2)),
      vecmem::get_data(*(this->m_output_host)));

  execute_sycl_test<matrix22_ops_functor<TypeParam>>(
      this->m_queue, this->m_m2->size(), vecmem::get_data(*(this->m_m2)),
      vecmem::get_data(*(this->m_output_device)));

  // Compare the outputs.
  this->compareOutputs();
}

/// Test for some operations with @c transform3
TYPED_TEST_P(sycl_transform_test, transform3D) {
  // Don't run the test at double precision, if the SYCL device doesn't
  // support it.
  if ((typeid(dscalar<TypeParam>) == typeid(double)) &&
      (this->m_queue.get_device().has(::sycl::aspect::fp64) == false)) {
    GTEST_SKIP();
  }

  // Run the test on the host, and on the/a device.
  execute_host_test<transform3_ops_functor<TypeParam>>(
      this->m_t1->size(), vecmem::get_data(*(this->m_t1)),
      vecmem::get_data(*(this->m_t2)), vecmem::get_data(*(this->m_t3)),
      vecmem::get_data(*(this->m_v1)), vecmem::get_data(*(this->m_v2)),
      vecmem::get_data(*(this->m_output_host)));

  execute_sycl_test<transform3_ops_functor<TypeParam>>(
      this->m_queue, this->m_t1->size(), vecmem::get_data(*(this->m_t1)),
      vecmem::get_data(*(this->m_t2)), vecmem::get_data(*(this->m_t3)),
      vecmem::get_data(*(this->m_v1)), vecmem::get_data(*(this->m_v2)),
      vecmem::get_data(*(this->m_output_device)));

  // Compare the outputs.
  this->compareOutputs();
}

}  // namespace detray::test::sycl
