// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/definitions/detail/qualifiers.hpp"

// Test include(s).
#include "detray/test/device/device_fixture.hpp"
#include "detray/test/framework/types.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/vector.hpp>
#include <vecmem/memory/memory_resource.hpp>

// GoogleTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <cmath>
#include <memory>

namespace detray::test {

/// Data fixture for matrix test cases
template <detray::concepts::algebra A>
class matrix_fixture : public device_fixture<dscalar<A>> {
  using index_t = dindex_type<A>;
  using scalar_t = dscalar<A>;
  template <index_t R, index_t C>
  using matrix_t = dmatrix<A, R, C>;

  using result_t = scalar_t;
  using base_fixture = device_fixture<result_t>;

 public:
  /// Constructor, setting up the inputs for all of the tests
  matrix_fixture(vecmem::memory_resource& mr) : base_fixture(mr) {}

 protected:
  /// Function setting things up for a test
  virtual void SetUp() override {
    // Set up the result vectors
    base_fixture::SetUp();

    // Set up the input vectors.
    m_m1 = std::make_unique<vecmem::vector<matrix_t<6, 4>>>(this->size(),
                                                            &this->resource());
    m_m2 = std::make_unique<vecmem::vector<matrix_t<2, 2>>>(this->size(),
                                                            &this->resource());

    // Initialise the input and output vectors.
    for (std::size_t i = 0u; i < this->size(); ++i) {
      for (std::size_t j = 0u; j < 6u; ++j) {
        for (std::size_t k = 0u; k < 4u; ++k) {
          detray::getter::element(m_m1->at(i), j, k) =
              static_cast<scalar_t>(j * 20.3 + k * 10.5);
        }
      }

      detray::getter::element(m_m2->at(i), 0u, 0u) = 4;
      detray::getter::element(m_m2->at(i), 0u, 1u) = 3;
      detray::getter::element(m_m2->at(i), 1u, 0u) = 12;
      detray::getter::element(m_m2->at(i), 1u, 1u) = 13;
    }
  }

  /// Function tearing things down after the test
  virtual void TearDown() override {
    // Delete the matrices.
    m_m1.reset();
    m_m2.reset();

    // Tear down the base class.
    base_fixture::TearDown();
  }

  /// @name Inputs for the tests
  /// @{
  std::unique_ptr<vecmem::vector<matrix_t<6, 4>>> m_m1;
  std::unique_ptr<vecmem::vector<matrix_t<2, 2>>> m_m2;
  /// @}

};  // class matrix_test_fixture

/// Functor running @c test_device_basics::matrix64_ops
template <detray::concepts::algebra A>
class matrix64_ops_functor {
  using index_t = dindex_type<A>;
  using scalar_t = dscalar<A>;
  template <index_t R, index_t C>
  using matrix_t = dmatrix<A, R, C>;

 public:
  DETRAY_HOST_DEVICE void operator()(
      std::size_t i, vecmem::data::vector_view<const matrix_t<6, 4>> m,
      vecmem::data::vector_view<scalar_t> output) const {
    // Create the VecMem vector(s).
    vecmem::device_vector<const matrix_t<6, 4>> vec_m(m);
    vecmem::device_vector<scalar_t> vec_output(output);

    // Perform the operation.
    auto ii = static_cast<typename decltype(vec_output)::size_type>(i);
    vec_output[ii] = matrix64_ops(vec_m[ii]);
  }

 private:
  /// Perform some trivial operations on an asymmetrix matrix
  DETRAY_HOST_DEVICE
  scalar_t matrix64_ops(const matrix_t<6, 4>& m) const {
    using namespace algebra;

    matrix_t<6, 4> m2;
    for (std::size_t i = 0u; i < 6u; ++i) {
      for (std::size_t j = 0u; j < 4u; ++j) {
        detray::getter::element(m2, i, j) = detray::getter::element(m, i, j);
      }
    }

    scalar_t result = 0.;
    for (std::size_t i = 0u; i < 6u; ++i) {
      for (std::size_t j = 0u; j < 4u; ++j) {
        result += 0.6f * detray::getter::element(m, i, j) +
                  0.7f * detray::getter::element(m2, i, j);
      }
    }

    // Test set_zero
    detray::matrix::set_zero(m2);
    for (std::size_t i = 0u; i < 6u; ++i) {
      for (std::size_t j = 0u; j < 4u; ++j) {
        result += 0.4f * detray::getter::element(m2, i, j);
      }
    }

    // Test set_identity
    detray::matrix::set_identity(m2);
    for (std::size_t i = 0u; i < 6u; ++i) {
      for (std::size_t j = 0u; j < 4u; ++j) {
        result += 0.3f * detray::getter::element(m2, i, j);
      }
    }

    // Test block operations
    auto b13 = detray::getter::block<1, 3>(m2, 0, 0);
    auto b13_tp = detray::matrix::transpose(b13);
    detray::getter::element(b13_tp, 0, 0) = 1;
    detray::getter::element(b13_tp, 1, 0) = 2;
    detray::getter::element(b13_tp, 2, 0) = 3;
    detray::getter::set_block(m2, b13_tp, 0, 0);

    auto b32 = detray::getter::block<3, 2>(m2, 2, 2);
    detray::getter::element(b32, 0, 0) = 4;
    detray::getter::element(b32, 0, 1) = 3;
    detray::getter::element(b32, 1, 0) = 12;
    detray::getter::element(b32, 1, 1) = 13;
    detray::getter::element(b32, 2, 0) = 5;
    detray::getter::element(b32, 2, 1) = 6;

    detray::getter::set_block(m2, b32, 2, 2);
    for (std::size_t i = 0u; i < 6u; ++i) {
      for (std::size_t j = 0u; j < 4u; ++j) {
        result += 0.57f * detray::getter::element(m2, i, j);
      }
    }

    return result;
  }
};

/// Functor running @c test_device_basics::matrix22_ops
template <detray::concepts::algebra A>
class matrix22_ops_functor {
  using index_t = dindex_type<A>;
  using scalar_t = dscalar<A>;
  template <index_t R, index_t C>
  using matrix_t = dmatrix<A, R, C>;

 public:
  DETRAY_HOST_DEVICE void operator()(
      std::size_t i, vecmem::data::vector_view<const matrix_t<2, 2>> m,
      vecmem::data::vector_view<scalar_t> output) const {
    // Create the VecMem vector(s).
    vecmem::device_vector<const matrix_t<2, 2>> vec_m(m);
    vecmem::device_vector<scalar_t> vec_output(output);

    // Perform the operation.
    auto ii = static_cast<typename decltype(vec_output)::size_type>(i);
    vec_output[ii] = matrix22_ops(vec_m[ii]);
  }

 private:
  /// Perform some trivial operations on an asymmetrix matrix
  DETRAY_HOST_DEVICE
  scalar_t matrix22_ops(const matrix_t<2, 2>& m22) const {
    // Test 2 X 2 matrix determinant
    auto m22_det = detray::matrix::determinant(m22);

    // Test 2 X 2 matrix inverse
    auto m22_inv = detray::matrix::inverse(m22);

    matrix_t<3, 3> m33;
    detray::getter::element(m33, 0, 0) = 1;
    detray::getter::element(m33, 0, 1) = 5;
    detray::getter::element(m33, 0, 2) = 7;
    detray::getter::element(m33, 1, 0) = 3;
    detray::getter::element(m33, 1, 1) = 5;
    detray::getter::element(m33, 1, 2) = 6;
    detray::getter::element(m33, 2, 0) = 2;
    detray::getter::element(m33, 2, 1) = 8;
    detray::getter::element(m33, 2, 2) = 9;

    // Test 3 X 3 matrix determinant
    auto m33_det = detray::matrix::determinant(m33);

    // Test 3 X 3 matrix inverse
    auto m33_inv = detray::matrix::inverse(m33);

    // Test Zero
    auto m23 = detray::matrix::zero<matrix_t<2, 3>>();
    detray::getter::element(m23, 0, 0) += 2;
    detray::getter::element(m23, 0, 1) += 3;
    detray::getter::element(m23, 0, 2) += 4;
    detray::getter::element(m23, 1, 0) += 5;
    detray::getter::element(m23, 1, 1) += 6;
    detray::getter::element(m23, 1, 2) += 7;

    // Test scalar X Matrix
    m23 = 2. * m23;

    // Test Transpose
    auto m32 = detray::matrix::transpose(m23);

    // Test Identity and (Matrix + Matrix)
    m32 = m32 + detray::matrix::identity<matrix_t<3, 2>>();

    // Test Matrix X scalar
    m32 = m32 * 2.;

    // Test Matrix multiplication
    auto new_m22 = m22_inv * m23 * m33_inv * m32;

    scalar_t result = 0;
    result += m22_det;
    result += m33_det;
    result += detray::getter::element(new_m22, 0, 0);
    result += detray::getter::element(new_m22, 0, 1);
    result += detray::getter::element(new_m22, 1, 0);
    result += detray::getter::element(new_m22, 1, 1);

    return result;
  }
};

}  // namespace detray::test
