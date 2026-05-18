// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Test include(s).
#include "detray/test/cpu/algebra_fixture.hpp"

// Project include(s)
#include "detray/algebra/utils/approximately_equal.hpp"
#include "detray/algebra/utils/casts.hpp"
#include "detray/algebra/utils/print.hpp"

// GoogleTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>

namespace detray::test {

#if !DETRAY_ALGEBRA_VC_AOS
TEST_F(detray_algebra, matrix_2x3) {
  static constexpr dindex_type<algebra_t> ROWS = 2;
  static constexpr dindex_type<algebra_t> COLS = 3;

  using matrix_2x3_t = dmatrix<algebra_t, 2, 3>;

  // Test type traits
  static_assert(std::is_same_v<detray::traits::index_t<matrix_2x3_t>,
                               dindex_type<algebra_t>>);
  static_assert(std::is_same_v<detray::traits::value_t<matrix_2x3_t>,
                               dscalar<algebra_t>>);
  static_assert(std::is_same_v<detray::traits::scalar_t<matrix_2x3_t>,
                               dscalar<algebra_t>>);
  static_assert(std::is_same_v<detray::traits::vector_t<matrix_2x3_t>,
                               detray::get_vector2D_t<algebra_t>>);

  static_assert(detray::traits::rows<matrix_2x3_t> == 2);
  static_assert(detray::traits::columns<matrix_2x3_t> == 3);
  static_assert(detray::traits::max_rank<matrix_2x3_t> == 2);
  static_assert(detray::traits::size<matrix_2x3_t> == 6);
  static_assert(!detray::traits::is_square<matrix_2x3_t>);
  static_assert(detray::traits::is_square<dmatrix<algebra_t, 2, 2>>);
  static_assert(detray::traits::is_square<dmatrix<algebra_t, 3, 3>>);

  // Test concepts
  static_assert(detray::concepts::matrix<matrix_2x3_t>);
  static_assert(!detray::concepts::scalar<matrix_2x3_t>);
  static_assert(!detray::concepts::vector<matrix_2x3_t>);
  static_assert(!detray::concepts::square_matrix<matrix_2x3_t>);

  static_assert(detray::concepts::index<detray::traits::index_t<matrix_2x3_t>>);
  static_assert(detray::concepts::value<detray::traits::value_t<matrix_2x3_t>>);
  static_assert(
      detray::concepts::scalar<detray::traits::scalar_t<matrix_2x3_t>>);
  static_assert(
      detray::concepts::vector<detray::traits::vector_t<matrix_2x3_t>>);

  // Test on matrix - vector operations
  detray::get_vector3D_t<algebra_t> vE{1.f, 2.f, 3.f};

  matrix_2x3_t m23;

  detray::getter::element(m23, 0, 0) = 1.f;
  detray::getter::element(m23, 0, 1) = 2.f;
  detray::getter::element(m23, 0, 2) = 3.f;
  detray::getter::element(m23, 1, 0) = 4.f;
  detray::getter::element(m23, 1, 1) = 5.f;
  detray::getter::element(m23, 1, 2) = 6.f;

  // Cast to (different) precision
  const auto m23_cast_f = detray::algebra::cast_to<float>(m23);

  for (std::size_t j = 0; j < 3; ++j) {
    for (std::size_t i = 0; i < 2; ++i) {
      auto elem_i = detray::getter::element(m23_cast_f, i, j);

      static_assert(std::same_as<decltype(elem_i), float>);
      ASSERT_FLOAT_EQ(elem_i,
                      static_cast<float>(detray::getter::element(m23, i, j)));
    }
  }

  const auto m23_cast_d = detray::algebra::cast_to<double>(m23);

  for (std::size_t j = 0; j < 3; ++j) {
    for (std::size_t i = 0; i < 2; ++i) {
      auto elem_i = detray::getter::element(m23_cast_d, i, j);

      static_assert(std::same_as<decltype(elem_i), double>);
      ASSERT_DOUBLE_EQ(elem_i,
                       static_cast<double>(detray::getter::element(m23, i, j)));
    }
  }

  const auto m23_cast_i = detray::algebra::cast_to<int>(m23);

  for (std::size_t j = 0; j < 3; ++j) {
    for (std::size_t i = 0; i < 2; ++i) {
      auto elem_i = detray::getter::element(m23_cast_i, i, j);

      static_assert(std::same_as<decltype(elem_i), int>);
      ASSERT_EQ(elem_i, static_cast<int>(detray::getter::element(m23, i, j)));
    }
  }

  detray::get_vector2D_t<algebra_t> v2 = m23 * vE;

  ASSERT_NEAR(v2[0], 14, this->epsilon());
  ASSERT_NEAR(v2[1], 32, this->epsilon());

  this->template matrix_test_ops_any_matrix<ROWS, COLS>();
}

TEST_F(detray_algebra, matrix_3x1) {
  static constexpr dindex_type<algebra_t> ROWS = 3;
  static constexpr dindex_type<algebra_t> COLS = 1;

  // Cross product on vector3 and matrix<3,1>
  dmatrix<algebra_t, 3, 1> vF;
  detray::getter::element(vF, 0, 0) = 5.f;
  detray::getter::element(vF, 1, 0) = 6.f;
  detray::getter::element(vF, 2, 0) = 13.f;

  // Test printing
  std::cout << vF << std::endl;

  // Cast to (different) precision
  const auto vF_cast_f = detray::algebra::cast_to<float>(vF);

  for (std::size_t i = 0; i < 3; ++i) {
    auto elem_i = detray::getter::element(vF_cast_f, i, 0);

    static_assert(std::same_as<decltype(elem_i), float>);
    ASSERT_FLOAT_EQ(elem_i,
                    static_cast<float>(detray::getter::element(vF, i, 0)));
  }

  const auto vF_cast_d = detray::algebra::cast_to<double>(vF);

  for (std::size_t i = 0; i < 3; ++i) {
    auto elem_i = detray::getter::element(vF_cast_d, i, 0);

    static_assert(std::same_as<decltype(elem_i), double>);
    ASSERT_DOUBLE_EQ(elem_i,
                     static_cast<double>(detray::getter::element(vF, i, 0)));
  }

  const auto vF_cast_i = detray::algebra::cast_to<int>(vF);

  for (std::size_t i = 0; i < 3; ++i) {
    auto elem_i = detray::getter::element(vF_cast_i, i, 0);

    static_assert(std::same_as<decltype(elem_i), int>);
    ASSERT_EQ(elem_i, static_cast<int>(detray::getter::element(vF, i, 0)));
  }

  detray::get_vector3D_t<algebra_t> vD{1.f, 1.f, 1.f};
  detray::get_vector3D_t<algebra_t> vG = vector::cross(vD, vF);
  ASSERT_NEAR(vG[0], 7.f, this->epsilon());
  ASSERT_NEAR(vG[1], -8.f, this->epsilon());
  ASSERT_NEAR(vG[2], 1.f, this->epsilon());

  // Dot product on vector3 and matrix<3,1>
  auto dot = vector::dot(vG, vF);
  ASSERT_NEAR(dot, 0.f, this->epsilon());

  this->template matrix_test_ops_any_matrix<ROWS, COLS>();
}

TEST_F(detray_algebra, matrix_6x4) {
  // Create the matrix.
  static constexpr dindex_type<algebra_t> ROWS = 6;
  static constexpr dindex_type<algebra_t> COLS = 4;
  using matrix_6x4_t = dmatrix<algebra_t, ROWS, COLS>;
  matrix_6x4_t m;

  // Test type traits
  static_assert(std::is_same_v<detray::traits::index_t<matrix_6x4_t>,
                               dindex_type<algebra_t>>);
  static_assert(std::is_same_v<detray::traits::value_t<matrix_6x4_t>,
                               dscalar<algebra_t>>);
  static_assert(std::is_same_v<detray::traits::scalar_t<matrix_6x4_t>,
                               dscalar<algebra_t>>);

  static_assert(detray::traits::rows<matrix_6x4_t> == 6);
  static_assert(detray::traits::columns<matrix_6x4_t> == 4);
  static_assert(detray::traits::max_rank<matrix_6x4_t> == 4);
  static_assert(detray::traits::size<matrix_6x4_t> == 24);
  static_assert(!detray::traits::is_square<matrix_6x4_t>);
  static_assert(detray::traits::is_square<dmatrix<algebra_t, 4, 4>>);
  static_assert(detray::traits::is_square<dmatrix<algebra_t, 6, 6>>);

  // Test concepts
  static_assert(detray::concepts::matrix<matrix_6x4_t>);
  static_assert(!detray::concepts::scalar<matrix_6x4_t>);
  static_assert(!detray::concepts::vector<matrix_6x4_t>);
  static_assert(!detray::concepts::square_matrix<matrix_6x4_t>);

  // Fill it.
  for (std::size_t i = 0; i < ROWS; ++i) {
    for (std::size_t j = 0; j < COLS; ++j) {
      detray::getter::element(m, i, j) =
          0.5f * static_cast<dscalar<algebra_t>>(i + j);
    }
  }

  // Check its content.
  const dmatrix<algebra_t, ROWS, COLS>& m_const_ref = m;
  for (std::size_t i = 0; i < ROWS; ++i) {
    for (std::size_t j = 0; j < COLS; ++j) {
      const dscalar<algebra_t> ref =
          0.5f * static_cast<dscalar<algebra_t>>(i + j);
      ASSERT_NEAR(detray::getter::element(m, i, j), ref, this->epsilon());
      ASSERT_NEAR(detray::getter::element(m_const_ref, i, j), ref,
                  this->epsilon());
    }
  }

  // Test set_zero
  detray::matrix::set_zero(m);
  for (std::size_t i = 0; i < ROWS; ++i) {
    for (std::size_t j = 0; j < COLS; ++j) {
      ASSERT_NEAR(detray::getter::element(m, i, j), 0., this->epsilon());
    }
  }

  // Test set_identity
  detray::matrix::set_identity(m);
  for (std::size_t i = 0; i < ROWS; ++i) {
    for (std::size_t j = 0; j < COLS; ++j) {
      if (i == j) {
        ASSERT_NEAR(detray::getter::element(m, i, j), 1.f, this->epsilon());
      } else {
        ASSERT_NEAR(detray::getter::element(m, i, j), 0.f, this->epsilon());
      }
    }
  }

  // Test block operations
  auto b13 = detray::getter::block<1, 3>(m, 0, 0);
  ASSERT_NEAR(detray::getter::element(b13, 0, 0), 1.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(b13, 0, 1), 0.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(b13, 0, 2), 0.f, this->epsilon());

  auto b13_tp = detray::matrix::transpose(b13);
  ASSERT_NEAR(detray::getter::element(b13_tp, 0, 0), 1.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(b13_tp, 1, 0), 0.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(b13_tp, 2, 0), 0.f, this->epsilon());

  auto b32 = detray::getter::block<3, 2>(m, 2, 2);
  ASSERT_NEAR(detray::getter::element(b32, 0, 0), 1.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(b32, 0, 1), 0.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(b32, 1, 0), 0.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(b32, 1, 1), 1.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(b32, 2, 0), 0.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(b32, 2, 1), 0.f, this->epsilon());

  detray::getter::element(b32, 0, 0) = 4.f;
  detray::getter::element(b32, 0, 1) = 3.f;
  detray::getter::element(b32, 1, 0) = 12.f;
  detray::getter::element(b32, 1, 1) = 13.f;
  detray::getter::element(b32, 2, 0) = 5.f;
  detray::getter::element(b32, 2, 1) = 6.f;

  detray::getter::set_block(m, b32, 2, 2);
  ASSERT_NEAR(detray::getter::element(m, 2, 2), 4.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(m, 2, 3), 3.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(m, 3, 2), 12.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(m, 3, 3), 13.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(m, 4, 2), 5.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(m, 4, 3), 6.f, this->epsilon());

  detray::get_vector3D_t<algebra_t> v = {10.f, 20.f, 30.f};
  detray::getter::set_block(m, v, 0, 2);
  ASSERT_NEAR(detray::getter::element(m, 0, 2), 10., this->epsilon());
  ASSERT_NEAR(detray::getter::element(m, 1, 2), 20., this->epsilon());
  ASSERT_NEAR(detray::getter::element(m, 2, 2), 30., this->epsilon());

  // Test printing
  std::cout << m << std::endl;

  this->template matrix_test_ops_any_matrix<ROWS, COLS>();
}

TEST_F(detray_algebra, matrix_3x3) {
  static constexpr dindex_type<algebra_t> N = 3;

  {
    detray::get_vector3D_t<algebra_t> v = {10.f, 20.f, 30.f};
    dmatrix<algebra_t, 3, 3> m33;
    detray::getter::element(m33, 0, 0) = 1;
    detray::getter::element(m33, 1, 0) = 2;
    detray::getter::element(m33, 2, 0) = 3;
    detray::getter::element(m33, 0, 1) = 5;
    detray::getter::element(m33, 1, 1) = 6;
    detray::getter::element(m33, 2, 1) = 7;
    detray::getter::element(m33, 0, 2) = 9;
    detray::getter::element(m33, 1, 2) = 10;
    detray::getter::element(m33, 2, 2) = 11;

    const detray::get_vector3D_t<algebra_t> v2 = m33 * v;
    ASSERT_NEAR(v2[0], 380., this->epsilon());
    ASSERT_NEAR(v2[1], 440., this->epsilon());
    ASSERT_NEAR(v2[2], 500., this->epsilon());
  }

  {
    dmatrix<algebra_t, 3, 3> m33;
    detray::getter::element(m33, 0, 0) = 1.f;
    detray::getter::element(m33, 0, 1) = 5.f;
    detray::getter::element(m33, 0, 2) = 7.f;
    detray::getter::element(m33, 1, 0) = 3.f;
    detray::getter::element(m33, 1, 1) = 5.f;
    detray::getter::element(m33, 1, 2) = 6.f;
    detray::getter::element(m33, 2, 0) = 2.f;
    detray::getter::element(m33, 2, 1) = 8.f;
    detray::getter::element(m33, 2, 2) = 9.f;

    // Test 3 X 3 matrix determinant
    auto m33_det = detray::matrix::determinant(m33);
    ASSERT_NEAR(m33_det, 20.f, this->tolerance());

    // Test 3 X 3 matrix inverse
    auto m33_inv = detray::matrix::inverse(m33);
    ASSERT_NEAR(detray::getter::element(m33_inv, 0, 0), -3.f / 20.f,
                this->tolerance());
    ASSERT_NEAR(detray::getter::element(m33_inv, 0, 1), 11.f / 20.f,
                this->tolerance());
    ASSERT_NEAR(detray::getter::element(m33_inv, 0, 2), -5.f / 20.f,
                this->tolerance());
    ASSERT_NEAR(detray::getter::element(m33_inv, 1, 0), -15.f / 20.f,
                this->tolerance());
    ASSERT_NEAR(detray::getter::element(m33_inv, 1, 1), -5.f / 20.f,
                this->tolerance());
    ASSERT_NEAR(detray::getter::element(m33_inv, 1, 2), 15.f / 20.f,
                this->tolerance());
    ASSERT_NEAR(detray::getter::element(m33_inv, 2, 0), 14.f / 20.f,
                this->tolerance());
    ASSERT_NEAR(detray::getter::element(m33_inv, 2, 1), 2.f / 20.f,
                this->tolerance());
    ASSERT_NEAR(detray::getter::element(m33_inv, 2, 2), -10.f / 20.f,
                this->tolerance());
  }

  this->template matrix_test_ops_square_matrix<N>();
}

TEST_F(detray_algebra, matrix_2x2) {
  static constexpr dindex_type<algebra_t> N = 2;

  dmatrix<algebra_t, 2, 2> m22;
  detray::getter::element(m22, 0, 0) = 4.f;
  detray::getter::element(m22, 0, 1) = 3.f;
  detray::getter::element(m22, 1, 0) = 12.f;
  detray::getter::element(m22, 1, 1) = 13.f;

  // Test 2 X 2 matrix determinant
  auto m22_det = detray::matrix::determinant(m22);
  ASSERT_NEAR(m22_det, 16.f, this->tolerance());

  // Test 2 X 2 matrix inverse
  auto m22_inv = detray::matrix::inverse(m22);
  ASSERT_NEAR(detray::getter::element(m22_inv, 0, 0), 13.f / 16.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m22_inv, 0, 1), -3.f / 16.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m22_inv, 1, 0), -12.f / 16.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m22_inv, 1, 1), 4.f / 16.f,
              this->tolerance());

  this->template matrix_test_ops_square_matrix<N>();
}

TEST_F(detray_algebra, matrix_5x5) {
  // Test 5 X 5 matrix
  dmatrix<algebra_t, 5, 5> m55;
  detray::getter::element(m55, 0, 0) = 1.f;
  detray::getter::element(m55, 0, 1) = 3.f;
  detray::getter::element(m55, 0, 2) = -9.f;
  detray::getter::element(m55, 0, 3) = -5.f;
  detray::getter::element(m55, 0, 4) = -2.f;

  detray::getter::element(m55, 1, 0) = -6.f;
  detray::getter::element(m55, 1, 1) = -3.f;
  detray::getter::element(m55, 1, 2) = 1.f;
  detray::getter::element(m55, 1, 3) = 0.f;
  detray::getter::element(m55, 1, 4) = 2.f;

  detray::getter::element(m55, 2, 0) = 12.f;
  detray::getter::element(m55, 2, 1) = 7.f;
  detray::getter::element(m55, 2, 2) = -9.f;
  detray::getter::element(m55, 2, 3) = 11.f;
  detray::getter::element(m55, 2, 4) = 2.f;

  detray::getter::element(m55, 3, 0) = -3.f;
  detray::getter::element(m55, 3, 1) = 4.f;
  detray::getter::element(m55, 3, 2) = 5.f;
  detray::getter::element(m55, 3, 3) = -6.f;
  detray::getter::element(m55, 3, 4) = 7.f;

  detray::getter::element(m55, 4, 0) = 9.f;
  detray::getter::element(m55, 4, 1) = 6.f;
  detray::getter::element(m55, 4, 2) = 3.f;
  detray::getter::element(m55, 4, 3) = 0.f;
  detray::getter::element(m55, 4, 4) = -3.f;

  auto m55_det = detray::matrix::determinant(m55);
  ASSERT_NEAR((m55_det - 17334.f) / 17334.f, 0.f, this->tolerance());

  auto m55_inv = detray::matrix::inverse(m55);

  ASSERT_NEAR(detray::getter::element(m55_inv, 0, 0), -2106.f / 17334.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m55_inv, 0, 1), -12312.f / 17334.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m55_inv, 0, 2), -486.f / 17334.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m55_inv, 0, 3), 864.f / 17334.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m55_inv, 0, 4), -5112.f / 17334.f,
              this->tolerance());

  ASSERT_NEAR(detray::getter::element(m55_inv, 1, 0), 2754.f / 17334.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m55_inv, 1, 1), 13878.f / 17334.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m55_inv, 1, 2), 1080.f / 17334.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m55_inv, 1, 3), -315.f / 17334.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m55_inv, 1, 4), 7401.f / 17334.f,
              this->tolerance());

  ASSERT_NEAR(detray::getter::element(m55_inv, 2, 0), -918.f / 17334.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m55_inv, 2, 1), 1152.f / 17334.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m55_inv, 2, 2), -360.f / 17334.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m55_inv, 2, 3), 105.f / 17334.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m55_inv, 2, 4), 1385.f / 17334.f,
              this->tolerance());

  ASSERT_NEAR(detray::getter::element(m55_inv, 3, 0), 108.f / 17334.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m55_inv, 3, 1), 7002.f / 17334.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m55_inv, 3, 2), 1062.f / 17334.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m55_inv, 3, 3), -1032.f / 17334.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m55_inv, 3, 4), 2896.f / 17334.f,
              this->tolerance());

  ASSERT_NEAR(detray::getter::element(m55_inv, 4, 0), -1728.f / 17334.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m55_inv, 4, 1), -8028.f / 17334.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m55_inv, 4, 2), 342.f / 17334.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m55_inv, 4, 3), 2067.f / 17334.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m55_inv, 4, 4), -4927.f / 17334.f,
              this->tolerance());
}

TEST_F(detray_algebra, matrix_6x6) {
  static constexpr dindex_type<algebra_t> N = 6;

  // Test 6 X 6 big matrix determinant
  dmatrix<algebra_t, 6, 6> m66_big;
  detray::getter::element(m66_big, 0, 0) = 1.f;
  detray::getter::element(m66_big, 0, 1) = 0.f;
  detray::getter::element(m66_big, 0, 2) = 3.f;
  detray::getter::element(m66_big, 0, 3) = 0.f;
  detray::getter::element(m66_big, 0, 4) = 0.f;
  detray::getter::element(m66_big, 0, 5) = 0.f;

  detray::getter::element(m66_big, 1, 0) = 0.f;
  detray::getter::element(m66_big, 1, 1) = -2.f;
  detray::getter::element(m66_big, 1, 2) = 4.f;
  detray::getter::element(m66_big, 1, 3) = 0.f;
  detray::getter::element(m66_big, 1, 4) = 5.f;
  detray::getter::element(m66_big, 1, 5) = 0.f;

  detray::getter::element(m66_big, 2, 0) = 0.f;
  detray::getter::element(m66_big, 2, 1) = 0.f;
  detray::getter::element(m66_big, 2, 2) = 3.f;
  detray::getter::element(m66_big, 2, 3) = 0.f;
  detray::getter::element(m66_big, 2, 4) = 0.f;
  detray::getter::element(m66_big, 2, 5) = 0.f;

  detray::getter::element(m66_big, 3, 0) = 0.f;
  detray::getter::element(m66_big, 3, 1) = 0.f;
  detray::getter::element(m66_big, 3, 2) = 0.f;
  detray::getter::element(m66_big, 3, 3) = 4.f;
  detray::getter::element(m66_big, 3, 4) = 0.f;
  detray::getter::element(m66_big, 3, 5) = 0.f;

  detray::getter::element(m66_big, 4, 0) = 0.f;
  detray::getter::element(m66_big, 4, 1) = 0.f;
  detray::getter::element(m66_big, 4, 2) = 0.f;
  detray::getter::element(m66_big, 4, 3) = 0.f;
  detray::getter::element(m66_big, 4, 4) = 9.f;
  detray::getter::element(m66_big, 4, 5) = 0.f;

  detray::getter::element(m66_big, 5, 0) = -1.f;
  detray::getter::element(m66_big, 5, 1) = -1.f;
  detray::getter::element(m66_big, 5, 2) = -1.f;
  detray::getter::element(m66_big, 5, 3) = -1.f;
  detray::getter::element(m66_big, 5, 4) = -1.f;
  detray::getter::element(m66_big, 5, 5) = -1.f;

  auto m66_big_det = detray::matrix::determinant(m66_big);
  ASSERT_NEAR(m66_big_det, 216.f, 2.f * this->tolerance());

  auto m66_big_inv = detray::matrix::inverse(m66_big);

  // Test comparison
  constexpr auto epsilon{std::numeric_limits<dscalar<algebra_t>>::epsilon()};

  ASSERT_TRUE(detray::algebra::approx_equal(m66_big, m66_big));
  ASSERT_TRUE(detray::algebra::approx_equal(m66_big, m66_big, epsilon));
  ASSERT_TRUE(detray::algebra::approx_equal(m66_big_inv, m66_big_inv));
  ASSERT_TRUE(detray::algebra::approx_equal(m66_big_inv, m66_big_inv, epsilon));

  dscalar<algebra_t> rel_err{1.f + 10.f * epsilon};
  dmatrix<algebra_t, 6, 6> m66_big_err = rel_err * m66_big;
  ASSERT_TRUE(
      detray::algebra::approx_equal(m66_big, m66_big_err, 11.f * epsilon));
  ASSERT_FALSE(
      detray::algebra::approx_equal(m66_big, m66_big_err, 9.f * epsilon));

  rel_err = 1.f + 17.f * epsilon;
  m66_big_err = rel_err * m66_big;
  ASSERT_TRUE(
      detray::algebra::approx_equal(m66_big, m66_big_err, 18.f * epsilon));
  ASSERT_FALSE(
      detray::algebra::approx_equal(m66_big, m66_big_err, 16.f * epsilon));

  dmatrix<algebra_t, 6, 6> prod66 = m66_big_inv * m66_big;
  auto I66 = detray::matrix::identity<decltype(m66_big)>();
  // Set the max error for Fastor and SMatrix plugins
  ASSERT_TRUE(
      detray::algebra::approx_equal(prod66, I66, epsilon, 2.f * epsilon));

  ASSERT_NEAR(detray::getter::element(m66_big_inv, 0, 0), 36.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 0, 1), 0.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 0, 2), -36.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 0, 3), 0.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 0, 4), 0.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 0, 5), 0.f / 36.f,
              this->tolerance());

  ASSERT_NEAR(detray::getter::element(m66_big_inv, 1, 0), 0.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 1, 1), -18.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 1, 2), 24.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 1, 3), 0.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 1, 4), 10.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 1, 5), 0.f / 36.f,
              this->tolerance());

  ASSERT_NEAR(detray::getter::element(m66_big_inv, 2, 0), 0.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 2, 1), 0.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 2, 2), 12.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 2, 3), 0.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 2, 4), 0.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 2, 5), 0.f / 36.f,
              this->tolerance());

  ASSERT_NEAR(detray::getter::element(m66_big_inv, 3, 0), 0.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 3, 1), 0.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 3, 2), 0.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 3, 3), 9.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 3, 4), 0.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 3, 5), 0.f / 36.f,
              this->tolerance());

  ASSERT_NEAR(detray::getter::element(m66_big_inv, 4, 0), 0.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 4, 1), 0.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 4, 2), 0.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 4, 3), 0.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 4, 4), 4.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 4, 5), 0.f / 36.f,
              this->tolerance());

  ASSERT_NEAR(detray::getter::element(m66_big_inv, 5, 0), -36.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 5, 1), 18.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 5, 2), 0.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 5, 3), -9.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 5, 4), -14.f / 36.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m66_big_inv, 5, 5), -36.f / 36.f,
              this->tolerance());

  // Cast to (different) precision
  const auto m66_cast_f = detray::algebra::cast_to<float>(m66_big_inv);

  for (std::size_t j = 0; j < 6; ++j) {
    for (std::size_t i = 0; i < 6; ++i) {
      auto elem_i = detray::getter::element(m66_cast_f, i, j);

      static_assert(std::same_as<decltype(elem_i), float>);
      ASSERT_FLOAT_EQ(elem_i, static_cast<float>(
                                  detray::getter::element(m66_big_inv, i, j)));
    }
  }

  const auto m66_cast_d = detray::algebra::cast_to<double>(m66_big_inv);

  for (std::size_t j = 0; j < 6; ++j) {
    for (std::size_t i = 0; i < 6; ++i) {
      auto elem_i = detray::getter::element(m66_cast_d, i, j);

      static_assert(std::same_as<decltype(elem_i), double>);
      ASSERT_DOUBLE_EQ(elem_i, static_cast<double>(
                                   detray::getter::element(m66_big_inv, i, j)));
    }
  }

  const auto m66_cast_i = detray::algebra::cast_to<int>(m66_big_inv);

  for (std::size_t j = 0; j < 6; ++j) {
    for (std::size_t i = 0; i < 6; ++i) {
      auto elem_i = detray::getter::element(m66_cast_i, i, j);

      static_assert(std::same_as<decltype(elem_i), int>);
      ASSERT_EQ(elem_i,
                static_cast<int>(detray::getter::element(m66_big_inv, i, j)));
    }
  }

  // Test 6 X 6 small matrix determinant
  dmatrix<algebra_t, 6, 6> m66_small;

  detray::getter::element(m66_small, 0, 0) =
      static_cast<dscalar<algebra_t>>(10.792386);
  detray::getter::element(m66_small, 0, 1) =
      static_cast<dscalar<algebra_t>>(0.216181);
  detray::getter::element(m66_small, 0, 2) =
      static_cast<dscalar<algebra_t>>(0.057650);
  detray::getter::element(m66_small, 0, 3) =
      static_cast<dscalar<algebra_t>>(-0.002764);
  detray::getter::element(m66_small, 0, 4) =
      static_cast<dscalar<algebra_t>>(0.000001);
  detray::getter::element(m66_small, 0, 5) = static_cast<dscalar<algebra_t>>(0);

  detray::getter::element(m66_small, 1, 0) =
      static_cast<dscalar<algebra_t>>(43.909368);
  detray::getter::element(m66_small, 1, 1) =
      static_cast<dscalar<algebra_t>>(10.372997);
  detray::getter::element(m66_small, 1, 2) =
      static_cast<dscalar<algebra_t>>(0.231496);
  detray::getter::element(m66_small, 1, 3) =
      static_cast<dscalar<algebra_t>>(-0.065972);
  detray::getter::element(m66_small, 1, 4) =
      static_cast<dscalar<algebra_t>>(-0.000002);
  detray::getter::element(m66_small, 1, 5) = static_cast<dscalar<algebra_t>>(0);

  detray::getter::element(m66_small, 2, 0) =
      static_cast<dscalar<algebra_t>>(0.045474);
  detray::getter::element(m66_small, 2, 1) =
      static_cast<dscalar<algebra_t>>(-0.001730);
  detray::getter::element(m66_small, 2, 2) =
      static_cast<dscalar<algebra_t>>(0.000246);
  detray::getter::element(m66_small, 2, 3) =
      static_cast<dscalar<algebra_t>>(0.000004);
  detray::getter::element(m66_small, 2, 4) = static_cast<dscalar<algebra_t>>(0);
  detray::getter::element(m66_small, 2, 5) = static_cast<dscalar<algebra_t>>(0);

  detray::getter::element(m66_small, 3, 0) =
      static_cast<dscalar<algebra_t>>(-0.255134);
  detray::getter::element(m66_small, 3, 1) =
      static_cast<dscalar<algebra_t>>(-0.059846);
  detray::getter::element(m66_small, 3, 2) =
      static_cast<dscalar<algebra_t>>(-0.001345);
  detray::getter::element(m66_small, 3, 3) =
      static_cast<dscalar<algebra_t>>(0.000383);
  detray::getter::element(m66_small, 3, 4) = static_cast<dscalar<algebra_t>>(0);
  detray::getter::element(m66_small, 3, 5) = static_cast<dscalar<algebra_t>>(0);

  detray::getter::element(m66_small, 4, 0) =
      static_cast<dscalar<algebra_t>>(-0.001490);
  detray::getter::element(m66_small, 4, 1) =
      static_cast<dscalar<algebra_t>>(-0.000057);
  detray::getter::element(m66_small, 4, 2) =
      static_cast<dscalar<algebra_t>>(-0.000008);
  detray::getter::element(m66_small, 4, 3) =
      static_cast<dscalar<algebra_t>>(0.000001);
  detray::getter::element(m66_small, 4, 4) =
      static_cast<dscalar<algebra_t>>(0.000001);
  detray::getter::element(m66_small, 4, 5) = static_cast<dscalar<algebra_t>>(0);

  detray::getter::element(m66_small, 5, 0) = static_cast<dscalar<algebra_t>>(0);
  detray::getter::element(m66_small, 5, 1) = static_cast<dscalar<algebra_t>>(0);
  detray::getter::element(m66_small, 5, 2) = static_cast<dscalar<algebra_t>>(0);
  detray::getter::element(m66_small, 5, 3) = static_cast<dscalar<algebra_t>>(0);
  detray::getter::element(m66_small, 5, 4) = static_cast<dscalar<algebra_t>>(0);
  detray::getter::element(m66_small, 5, 5) =
      static_cast<dscalar<algebra_t>>(89875.517874);

  auto m66_small_det = detray::matrix::determinant(m66_small);
  ASSERT_NEAR((m66_small_det - 4.30636e-11f) / 4.30636e-11f, 0.f,
              2.f * this->tolerance());

  this->template matrix_test_ops_square_matrix<N>();
}

TEST_F(detray_algebra, matrix_7x4_4x12) {
  this->template matrix_test_ops_inhomogeneous_multipliable_matrices<7, 4,
                                                                     12>();
}

TEST_F(detray_algebra, matrix_17x9_9x4) {
  this->template matrix_test_ops_inhomogeneous_multipliable_matrices<17, 9,
                                                                     4>();
}

TEST_F(detray_algebra, matrix_5x2_2x3) {
  this->template matrix_test_ops_inhomogeneous_multipliable_matrices<5, 2, 3>();
}

TEST_F(detray_algebra, matrix_small_mixed) {
  dmatrix<algebra_t, 2, 2> m22;
  detray::getter::element(m22, 0, 0) = 4.f;
  detray::getter::element(m22, 0, 1) = 3.f;
  detray::getter::element(m22, 1, 0) = 12.f;
  detray::getter::element(m22, 1, 1) = 13.f;

  // Test 2 X 2 matrix inverse
  auto m22_inv = detray::matrix::inverse(m22);
  ASSERT_NEAR(detray::getter::element(m22_inv, 0, 0), 13.f / 16.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m22_inv, 0, 1), -3.f / 16.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m22_inv, 1, 0), -12.f / 16.f,
              this->tolerance());
  ASSERT_NEAR(detray::getter::element(m22_inv, 1, 1), 4.f / 16.f,
              this->tolerance());

  dmatrix<algebra_t, 3, 3> m33;
  detray::getter::element(m33, 0, 0) = 1.f;
  detray::getter::element(m33, 0, 1) = 5.f;
  detray::getter::element(m33, 0, 2) = 7.f;
  detray::getter::element(m33, 1, 0) = 3.f;
  detray::getter::element(m33, 1, 1) = 5.f;
  detray::getter::element(m33, 1, 2) = 6.f;
  detray::getter::element(m33, 2, 0) = 2.f;
  detray::getter::element(m33, 2, 1) = 8.f;
  detray::getter::element(m33, 2, 2) = 9.f;

  auto m33_inv = detray::matrix::inverse(m33);

  // Test Zero
  auto m23 = detray::matrix::zero<dmatrix<algebra_t, 2, 3>>();
  detray::getter::element(m23, 0, 0) += 2.f;
  detray::getter::element(m23, 0, 1) += 3.f;
  detray::getter::element(m23, 0, 2) += 4.f;
  detray::getter::element(m23, 1, 0) += 5.f;
  detray::getter::element(m23, 1, 1) += 6.f;
  detray::getter::element(m23, 1, 2) += 7.f;

  // Test scalar X Matrix
  m23 = 2. * m23;
  ASSERT_NEAR(detray::getter::element(m23, 0, 0), 4.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(m23, 0, 1), 6.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(m23, 0, 2), 8.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(m23, 1, 0), 10.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(m23, 1, 1), 12.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(m23, 1, 2), 14.f, this->epsilon());

  // Test Transpose
  auto m32 = detray::matrix::transpose(m23);

  // Test Identity and (Matrix + Matrix)
  m32 = m32 + detray::matrix::identity<dmatrix<algebra_t, 3, 2>>();

  // Test Matrix X scalar
  m32 = m32 * 2.f;

  ASSERT_NEAR(detray::getter::element(m32, 0, 0), 10.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(m32, 0, 1), 20.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(m32, 1, 0), 12.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(m32, 1, 1), 26.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(m32, 2, 0), 16.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(m32, 2, 1), 28.f, this->epsilon());

  // Test Matrix multiplication
  m22 = m22_inv * m23 * m33_inv * m32;

  ASSERT_NEAR(detray::getter::element(m22, 0, 0), 6.225f, this->tolerance());
  ASSERT_NEAR(detray::getter::element(m22, 0, 1), 14.675f, this->tolerance());
  ASSERT_NEAR(detray::getter::element(m22, 1, 0), -3.3f, this->tolerance());
  ASSERT_NEAR(detray::getter::element(m22, 1, 1), -7.9f, this->tolerance());
}
#endif
}  // namespace detray::test
