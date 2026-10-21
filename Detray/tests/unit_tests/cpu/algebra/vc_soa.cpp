// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
// clang-format off
#include "detray/definitions/algebra.hpp"
#include "detray/algebra/utils/approximately_equal.hpp"
#include "detray/algebra/utils/casts.hpp"
#include "detray/algebra/utils/print.hpp"
// clang-format on

// Test include(s)
#include "detray/test/framework/types.hpp"

// GoogleTest include(s).
#include <gtest/gtest.h>

// System include(s)
#include <array>
#include <limits>

using namespace detray;

using value_t = float;

constexpr value_t tol{1e-5f};

/// This test the vector functions on an SoA (Vc::Vector) based vector
TEST(detray_algebra_soa, vector) {
  using test_algebra_t = detray::vc_soa<value_t>;
  using vector3_v = dvector3D<test_algebra_t>;

  // Value type is Vc::Vector<float>
  using scalar_t = dscalar<test_algebra_t>;

  static_assert(detray::concepts::scalar<scalar_t>);
  static_assert(detray::concepts::vector<vector3_v>);
  static_assert(detray::concepts::vector3D<vector3_v>);

  // Cast simd scalar to different precisions
  using scalar_f = Vc::float_v;
  using scalar_d = Vc::double_v;
  using scalar_i = Vc::int_v;

  auto s1 = detray::algebra::cast_to<float>(scalar_t(1.f));
  auto s2 = detray::algebra::cast_to<double>(scalar_t(2.f));
  auto s3 = detray::algebra::cast_to<int>(scalar_t(3.f));

  static_assert(std::same_as<decltype(s1), scalar_f>);
  static_assert(std::same_as<decltype(s2), scalar_d>);
  static_assert(std::same_as<decltype(s3), scalar_i>);

  ASSERT_TRUE((s1 == scalar_f(1.f)).isFull());
  ASSERT_TRUE((s2 == scalar_d(2.f)).isFull());
  ASSERT_TRUE((s3 == scalar_i(3.f)).isFull());

  vector3_v a{1.f, 2.f, 3.f};
  vector3_v b{4.f, 5.f, 6.f};

  EXPECT_TRUE((a[0] == scalar_t(1.f)).isFull());
  EXPECT_TRUE((a[1] == scalar_t(2.f)).isFull());
  EXPECT_TRUE((a[2] == scalar_t(3.f)).isFull());

  // Test printing
  std::cout << a << std::endl;

  // Test comparison
  constexpr auto epsilon{std::numeric_limits<value_t>::epsilon()};

  EXPECT_TRUE(detray::algebra::approx_equal(a, a));
  EXPECT_TRUE(detray::algebra::approx_equal(a, a, epsilon));
  EXPECT_FALSE(detray::algebra::approx_equal(a, b));

  value_t rel_err = 1.f + 10.f * epsilon;
  vector3_v a_err = rel_err * a;
  EXPECT_TRUE(detray::algebra::approx_equal(a, a_err, 11.f * epsilon));
  EXPECT_FALSE(detray::algebra::approx_equal(a, a_err, 9.f * epsilon));

  rel_err = 1.f + 17.f * epsilon;
  a_err = rel_err * a;
  EXPECT_TRUE(detray::algebra::approx_equal(a, a_err, 18.f * epsilon));
  EXPECT_FALSE(detray::algebra::approx_equal(a, a_err, 16.f * epsilon));

  // Swap an element
  vector3_v a_err_cpy = a_err;
  EXPECT_TRUE(a_err_cpy == a_err);
  EXPECT_TRUE(detray::algebra::approx_equal(a_err_cpy, a_err));

  auto& vec_elem = a_err[0];
  vec_elem[0] += 1.f;
  EXPECT_FALSE(a_err_cpy == a_err);
  EXPECT_FALSE(detray::algebra::approx_equal(a_err_cpy, a_err));
  // Cast simd vectors to different precision
  auto a_cast_f = detray::algebra::cast_to<float>(a);
  auto a_cast_d = detray::algebra::cast_to<double>(a);
  auto a_cast_i = detray::algebra::cast_to<int>(a);

  using algebra_f_t = detray::vc_soa<float>;
  using algebra_d_t = detray::vc_soa<double>;
  using algebra_i_t = detray::vc_soa<int>;
  static_assert(std::same_as<decltype(a_cast_f), dvector3D<algebra_f_t>>);
  static_assert(std::same_as<decltype(a_cast_d), dvector3D<algebra_d_t>>);
  static_assert(std::same_as<decltype(a_cast_i), dvector3D<algebra_i_t>>);

  for (int i = 0; i < 3; ++i) {
    EXPECT_TRUE(
        (a_cast_f[i] == detray::algebra::cast_to<float>(a[i])).isFull());
    EXPECT_TRUE(
        (a_cast_d[i] == detray::algebra::cast_to<double>(a[i])).isFull());
    EXPECT_TRUE((a_cast_i[i] == detray::algebra::cast_to<int>(a[i])).isFull());
  }

  // Masked comparison
  auto m = a.compare(a);
  EXPECT_TRUE(m[0].isFull());
  EXPECT_TRUE(m[1].isFull());
  EXPECT_TRUE(m[2].isFull());

  m = a.compare(b);
  EXPECT_FALSE(m[0].isFull());
  EXPECT_FALSE(m[1].isFull());
  EXPECT_FALSE(m[2].isFull());

  // Full comparisons
  EXPECT_TRUE(a == a);
  EXPECT_FALSE(a == b);

  // Addition
  auto v_add = a + b;
  EXPECT_TRUE((v_add[0] == scalar_t(5.f)).isFull());
  EXPECT_TRUE((v_add[1] == scalar_t(7.f)).isFull());
  EXPECT_TRUE((v_add[2] == scalar_t(9.f)).isFull());

  // Subration
  auto v_sub = a - b;
  EXPECT_TRUE((v_sub[0] == scalar_t(-3.f)).isFull());
  EXPECT_TRUE((v_sub[1] == scalar_t(-3.f)).isFull());
  EXPECT_TRUE((v_sub[2] == scalar_t(-3.f)).isFull());

  // Multiplication
  auto v_mul = a * b;
  EXPECT_TRUE((v_mul[0] == scalar_t(4.f)).isFull());
  EXPECT_TRUE((v_mul[1] == scalar_t(10.f)).isFull());
  EXPECT_TRUE((v_mul[2] == scalar_t(18.f)).isFull());

  // Division
  auto v_div = a / b;
  EXPECT_TRUE((v_div[0] == scalar_t(0.25f)).isFull());
  EXPECT_TRUE((v_div[1] == scalar_t(0.4f)).isFull());
  EXPECT_TRUE((v_div[2] == scalar_t(0.5f)).isFull());

  // Scalar multiplication
  auto v_smul = 2.f * b;
  EXPECT_TRUE((v_smul[0] == scalar_t(8.f)).isFull());
  EXPECT_TRUE((v_smul[1] == scalar_t(10.f)).isFull());
  EXPECT_TRUE((v_smul[2] == scalar_t(12.f)).isFull());

  // Expression
  auto v_expr = (b / a) - (2.5f * b) + vector3_v{};
  EXPECT_TRUE((v_expr[0] == scalar_t(-6.f)).isFull());
  EXPECT_TRUE((v_expr[1] == scalar_t(-10.f)).isFull());
  EXPECT_TRUE((v_expr[2] == scalar_t(-13.f)).isFull());

  auto d{vector::dot(a, b)};
  EXPECT_TRUE((d == scalar_t(32.f)).isFull());

  scalar_t norms_a{vector::norm(vector::normalize(a))};
  scalar_t norms_b{vector::norm(vector::normalize(b))};
  for (unsigned int i{0u}; i < norms_a.size(); ++i) {
    EXPECT_NEAR(norms_a[i], 1.f, tol);
    EXPECT_NEAR(norms_b[i], 1.f, tol);
  }

  auto cr{vector::cross(a, b)};
  EXPECT_TRUE((cr[0] == scalar_t(-3.f)).isFull());
  EXPECT_TRUE((cr[1] == scalar_t(6.f)).isFull());
  EXPECT_TRUE((cr[2] == scalar_t(-3.f)).isFull());

  static_assert(std::is_convertible_v<decltype(v_expr), vector3_v>,
                "expression type not convertible");
}

/// This test the getter functions on an SoA (Vc::Vector) based vector
TEST(detray_algebra_soa, getter) {
  using test_algebra_t = detray::vc_soa<value_t>;

  using vector3_v = dvector3D<test_algebra_t>;

  vector3_v a{1.f, 2.f, 3.f};

  // All results in the vector are the same, so only check the first one

  // Phi angle
  auto v_phi = vector::phi(a);
  EXPECT_NEAR(v_phi[0], static_cast<value_t>(std::atan2(2., 1.)), tol);

  // Perpendicular projection
  auto v_perp = vector::perp(a);
  EXPECT_NEAR(v_perp[0], std::sqrt(5.), tol);

  // Theta angle
  auto v_theta = vector::theta(a);
  EXPECT_NEAR(v_theta[0], static_cast<value_t>(std::atan2(std::sqrt(5.), 3.)),
              tol);

  // Norm of the vector
  auto v_norm = vector::norm(a);
  EXPECT_NEAR(v_norm[0], std::sqrt(14.), tol);

  // Eta of the vector
  auto v_eta = vector::eta(a);
  EXPECT_NEAR(v_eta[0],
              static_cast<value_t>(std::atanh(1. / std::sqrt(14.) * 3.)), tol);
}

/// This test an SoA (Vc::Vector) based affine transform3
TEST(detray_algebra_soa, transform3) {
  using test_algebra_t = detray::vc_soa<value_t>;

  // Print the linear algebra types of this backend
  using algebra::operator<<;

  using vector3 = dvector3D<test_algebra_t>;
  using point3 = dpoint3D<test_algebra_t>;
  using scalar_t = dscalar<test_algebra_t>;
  using transform3 = dtransform3D<test_algebra_t>;

  static_assert(detray::concepts::transform3D<transform3>);

  transform3 idty{};

  EXPECT_TRUE((idty(0, 0) == scalar_t::One()).isFull());
  EXPECT_TRUE((idty(1, 0) == scalar_t::Zero()).isFull());
  EXPECT_TRUE((idty(2, 0) == scalar_t::Zero()).isFull());
  EXPECT_TRUE((idty(0, 1) == scalar_t::Zero()).isFull());
  EXPECT_TRUE((idty(1, 1) == scalar_t::One()).isFull());
  EXPECT_TRUE((idty(2, 1) == scalar_t::Zero()).isFull());
  EXPECT_TRUE((idty(0, 2) == scalar_t::Zero()).isFull());
  EXPECT_TRUE((idty(1, 2) == scalar_t::Zero()).isFull());
  EXPECT_TRUE((idty(2, 2) == scalar_t::One()).isFull());
  EXPECT_TRUE((idty(0, 3) == scalar_t::Zero()).isFull());
  EXPECT_TRUE((idty(1, 3) == scalar_t::Zero()).isFull());
  EXPECT_TRUE((idty(2, 3) == scalar_t::Zero()).isFull());

  // Preparatioon work
  vector3 z = vector::normalize(vector3{3.f, 2.f, 1.f});
  vector3 x = vector::normalize(vector3{2.f, -3.f, 0.f});
  vector3 y = vector::cross(z, x);
  point3 t = {2.f, 3.f, 4.f};

  // Test constructor from t, z, x
  transform3 trf1(t, z, x);
  ASSERT_TRUE(trf1 == trf1);
  transform3 trf2;
  trf2 = trf1;

  // Test printing
  std::cout << trf1 << std::endl;

  // Test comparison
  constexpr auto epsilon{std::numeric_limits<value_t>::epsilon()};

  EXPECT_TRUE(detray::algebra::approx_equal(trf1, trf1));
  EXPECT_TRUE(detray::algebra::approx_equal(trf1, trf1, epsilon));

  value_t rel_err{1.f + 10.f * epsilon};
  transform3 trf1_err(rel_err * t, rel_err * z, rel_err * x);
  EXPECT_FALSE(trf1 == trf1_err);
  EXPECT_TRUE(detray::algebra::approx_equal(trf1, trf1_err, 200.f * epsilon));
  EXPECT_FALSE(detray::algebra::approx_equal(trf1, trf1_err, 10.f * epsilon));
  // Cast simd vectors to different precision
  auto trf1_cast_f = detray::algebra::cast_to<float>(trf1);
  auto trf1_cast_d = detray::algebra::cast_to<double>(trf1);
  auto trf1_cast_i = detray::algebra::cast_to<int>(trf1);

  using algebra_f_t = detray::vc_soa<float>;
  using algebra_d_t = detray::vc_soa<double>;
  using algebra_i_t = detray::vc_soa<int>;
  static_assert(std::same_as<decltype(trf1_cast_f), dtransform3D<algebra_f_t>>);
  static_assert(std::same_as<decltype(trf1_cast_d), dtransform3D<algebra_d_t>>);
  static_assert(std::same_as<decltype(trf1_cast_i), dtransform3D<algebra_i_t>>);

  const auto& mat_f = trf1_cast_f.matrix();
  const auto& mat_d = trf1_cast_d.matrix();
  const auto& mat_i = trf1_cast_i.matrix();
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      const auto& elem_ij = trf1.matrix()[i][j];
      EXPECT_TRUE(
          (mat_f[i][j] == detray::algebra::cast_to<float>(elem_ij)).isFull());
      EXPECT_TRUE(
          (mat_d[i][j] == detray::algebra::cast_to<double>(elem_ij)).isFull());
      EXPECT_TRUE(
          (mat_i[i][j] == detray::algebra::cast_to<int>(elem_ij)).isFull());
    }
  }

  EXPECT_TRUE((trf2(0, 0) == x[0]).isFull());
  EXPECT_TRUE((trf2(1, 0) == x[1]).isFull());
  EXPECT_TRUE((trf2(2, 0) == x[2]).isFull());
  EXPECT_TRUE((trf2(0, 1) == y[0]).isFull());
  EXPECT_TRUE((trf2(1, 1) == y[1]).isFull());
  EXPECT_TRUE((trf2(2, 1) == y[2]).isFull());
  EXPECT_TRUE((trf2(0, 2) == z[0]).isFull());
  EXPECT_TRUE((trf2(1, 2) == z[1]).isFull());
  EXPECT_TRUE((trf2(2, 2) == z[2]).isFull());
  EXPECT_TRUE((trf2(0, 3) == 2.f * scalar_t::One()).isFull());
  EXPECT_TRUE((trf2(1, 3) == 3.f * scalar_t::One()).isFull());
  EXPECT_TRUE((trf2(2, 3) == 4.f * scalar_t::One()).isFull());

  // Check that local origin translates into global translation
  point3 lzero = {0.f, 0.f, 0.f};
  point3 gzero = trf2.point_to_global(lzero);
  EXPECT_TRUE((gzero[0] == t[0]).isFull());
  EXPECT_TRUE((gzero[1] == t[1]).isFull());
  EXPECT_TRUE((gzero[2] == t[2]).isFull());

  // Check a round trip for point
  point3 loc_pt = {3.f, 4.f, 5.f};
  point3 glob_pt = trf2.point_to_global(loc_pt);
  point3 loc_pt_r = trf2.point_to_local(glob_pt);
  EXPECT_NEAR(loc_pt[0][0], loc_pt_r[0][0], tol);
  EXPECT_NEAR(loc_pt[1][0], loc_pt_r[1][0], tol);
  EXPECT_NEAR(loc_pt[2][0], loc_pt_r[2][0], tol);

  // Check a point versus vector transform
  // vector should not change if transformed by a pure translation
  transform3 ttrf(t);

  vector3 glob_vec = {1.f, 1.f, 1.f};
  vector3 loc_vec = ttrf.vector_to_local(glob_vec);
  EXPECT_NEAR(glob_vec[0][0], loc_vec[0][0], tol);
  EXPECT_NEAR(glob_vec[1][0], loc_vec[1][0], tol);
  EXPECT_NEAR(glob_vec[2][0], loc_vec[2][0], tol);

  // Check a round trip for vector
  vector3 loc_vecB = {7.f, 8.f, 9.f};
  vector3 glob_vecB = trf2.vector_to_local(loc_vecB);
  vector3 loc_vecC = trf2.vector_to_global(glob_vecB);
  EXPECT_NEAR(loc_vecB[0][0], loc_vecC[0][0], tol);
  EXPECT_NEAR(loc_vecB[1][0], loc_vecC[1][0], tol);
  EXPECT_NEAR(loc_vecB[2][0], loc_vecC[2][0], tol);
}

/// This test an SoA (Vc::Vector) based 2x3 matrix
TEST(detray_algebra_soa, matrix3) {
  using test_algebra_t = detray::vc_soa<value_t>;
  using matrix_2x3_t = dmatrix<test_algebra_t, 2, 3>;

  // Test type traits
  static_assert(
      std::is_same_v<detray::traits::index_t<matrix_2x3_t>, std::size_t>);
  static_assert(std::is_same_v<detray::traits::value_t<matrix_2x3_t>, value_t>);
  static_assert(std::is_same_v<detray::traits::scalar_t<matrix_2x3_t>,
                               Vc::Vector<value_t>>);
  static_assert(std::is_same_v<detray::traits::vector_t<matrix_2x3_t>,
                               dvector2D<test_algebra_t>>);

  static_assert(detray::traits::rows<matrix_2x3_t> == 2);
  static_assert(detray::traits::columns<matrix_2x3_t> == 3);
  static_assert(detray::traits::max_rank<matrix_2x3_t> == 2);
  static_assert(detray::traits::size<matrix_2x3_t> == 6);
  static_assert(!detray::traits::is_square<matrix_2x3_t>);
  static_assert(detray::traits::is_square<dmatrix<test_algebra_t, 2, 2>>);
  static_assert(detray::traits::is_square<dmatrix<test_algebra_t, 3, 3>>);
}

/// This test an SoA (Vc::Vector) based 6x4 matrix
TEST(detray_algebra_soa, matrix64) {
  using test_algebra_t = detray::vc_soa<value_t>;

  // Create the matrix.
  using matrix_6x4_t = dmatrix<test_algebra_t, 6, 4>;
  using scalar_t = dscalar<test_algebra_t>;

  matrix_6x4_t m;

  // Test type traits
  static_assert(
      std::is_same_v<detray::traits::index_t<matrix_6x4_t>, std::size_t>);
  static_assert(std::is_same_v<detray::traits::value_t<matrix_6x4_t>, value_t>);
  static_assert(std::is_same_v<detray::traits::scalar_t<matrix_6x4_t>,
                               Vc::Vector<value_t>>);

  static_assert(detray::traits::rows<matrix_6x4_t> == 6);
  static_assert(detray::traits::columns<matrix_6x4_t> == 4);
  static_assert(detray::traits::max_rank<matrix_6x4_t> == 4);
  static_assert(detray::traits::size<matrix_6x4_t> == 24);
  static_assert(!detray::traits::is_square<matrix_6x4_t>);
  static_assert(detray::traits::is_square<dmatrix<test_algebra_t, 4, 4>>);
  static_assert(detray::traits::is_square<dmatrix<test_algebra_t, 6, 6>>);

  // Test printing
  std::cout << m << std::endl;

  auto I64 = detray::matrix::identity<matrix_6x4_t>();

  // Test comparison
  constexpr auto epsilon{std::numeric_limits<value_t>::epsilon()};

  EXPECT_TRUE(detray::algebra::approx_equal(m, m));
  EXPECT_TRUE(detray::algebra::approx_equal(m, m, epsilon));
  EXPECT_FALSE(detray::algebra::approx_equal(m, I64));

  value_t rel_err{1.f + 10.f * epsilon};
  matrix_6x4_t I64_err = scalar_t(rel_err) * I64;
  EXPECT_FALSE(I64 == I64_err);
  EXPECT_TRUE(detray::algebra::approx_equal(I64, I64_err, 11.f * epsilon));
  EXPECT_FALSE(detray::algebra::approx_equal(I64, I64_err, 9.f * epsilon));
  // Cast simd vectors to different precision
  auto m_cast_f = detray::algebra::cast_to<float>(m);
  auto m_cast_d = detray::algebra::cast_to<double>(m);
  auto m_cast_i = detray::algebra::cast_to<int>(m);

  using algebra_f_t = detray::vc_soa<float>;
  using algebra_d_t = detray::vc_soa<double>;
  using algebra_i_t = detray::vc_soa<int>;
  static_assert(std::same_as<decltype(m_cast_f), dmatrix<algebra_f_t, 6, 4>>);
  static_assert(std::same_as<decltype(m_cast_d), dmatrix<algebra_d_t, 6, 4>>);
  static_assert(std::same_as<decltype(m_cast_i), dmatrix<algebra_i_t, 6, 4>>);
}
