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
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/rectangle2D.hpp"
#include "detray/geometry/surface_descriptor.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/utils/quadratic_equation.hpp"

// Test include(s)
#include "detray/test/framework/types.hpp"

// GoogleTest include(s).
#include <gtest/gtest.h>

// System include(s)
#include <array>
#include <cmath>
#include <limits>

using namespace detray;

using value_t = float;

constexpr value_t tol{1e-5f};

/// This tests the lane bundle that serves as the scalar of the SoA plugin
TEST(detray_algebra_soa, simd) {
  using test_algebra_t = detray::array_soa<value_t>;
  using scalar_t = dscalar<test_algebra_t>;
  using bool_t = dbool<test_algebra_t>;

  static_assert(detray::concepts::scalar<scalar_t>);
  static_assert(detray::concepts::simd_scalar<scalar_t>);
  static_assert(detray::concepts::soa<test_algebra_t>);
  static_assert(!detray::concepts::aos<test_algebra_t>);
  static_assert(detray::concepts::algebra<test_algebra_t>);
  static_assert(scalar_t::size() == 8u);
  static_assert(std::same_as<dscalar<detray::array_soa<value_t, 4>>,
                             algebra::array_soa::simd<value_t, 4>>);

  // Default construction is zero, like Vc::Vector
  scalar_t z;
  EXPECT_TRUE(detray::detail::all_of(z == detray::detail::zero<scalar_t>()));
  EXPECT_TRUE(detray::detail::all_of(z == 0.f));

  // Broadcast and lane index
  scalar_t two{2.f};
  const scalar_t idx = detray::detail::iota<scalar_t>();
  for (std::size_t i = 0u; i < scalar_t::size(); ++i) {
    EXPECT_EQ(two[i], 2.f);
    EXPECT_EQ(idx[i], static_cast<value_t>(i));
  }

  // Arithmetic with bundles and with single values on both sides
  const scalar_t a = idx + 1.f;
  const scalar_t b = 2.f * a;
  EXPECT_TRUE(detray::detail::all_of(b == a + a));
  EXPECT_TRUE(detray::detail::all_of(b - a == a));
  EXPECT_TRUE(detray::detail::all_of(b / 2.f == a));
  EXPECT_TRUE(detray::detail::all_of(-a == 0.f - a));
  EXPECT_TRUE(detray::detail::all_of(a * b == b * a));
  EXPECT_TRUE(detray::detail::all_of(1.f / detray::detail::one<scalar_t>() ==
                                     detray::detail::one<scalar_t>()));

  scalar_t c = a;
  c += 1.f;
  c -= detray::detail::one<scalar_t>();
  c *= 2.f;
  c /= two;
  EXPECT_TRUE(detray::detail::all_of(c == a));

  // Comparisons yield masks
  const bool_t m = idx < 3.f;
  static_assert(std::same_as<decltype(m), const bool_t>);
  EXPECT_FALSE(detray::detail::all_of(m));
  EXPECT_FALSE(detray::detail::none_of(m));
  EXPECT_TRUE(detray::detail::any_of(m));
  EXPECT_EQ(detray::detail::count(m), 3u);
  EXPECT_TRUE(detray::detail::count(!m) == scalar_t::size() - 3u);
  EXPECT_TRUE(detray::detail::none_of(m && !m));
  EXPECT_TRUE(detray::detail::all_of(m || !m));
  EXPECT_TRUE(detray::detail::all_of(idx >= 0.f));
  EXPECT_TRUE((3.f > idx) == m);
  EXPECT_TRUE(detray::detail::all_of(bool_t{true}));
  EXPECT_TRUE(detray::detail::none_of(bool_t{false}));

  // Boolean reductions used by the core
  EXPECT_TRUE(detray::detail::any_of(m));
  EXPECT_FALSE(detray::detail::all_of(m));
  EXPECT_FALSE(detray::detail::none_of(m));
  EXPECT_TRUE(detray::detail::none_of(idx < 0.f));

  // Masked assignment
  scalar_t d = idx;
  detray::detail::set_if(d, m, scalar_t{42.f});
  for (std::size_t i = 0u; i < scalar_t::size(); ++i) {
    EXPECT_EQ(d[i], i < 3u ? 42.f : static_cast<value_t>(i));
  }
  detray::detail::set_if(d, !m, a);
  for (std::size_t i = 0u; i < scalar_t::size(); ++i) {
    EXPECT_EQ(d[i], i < 3u ? 42.f : static_cast<value_t>(i + 1u));
  }

  scalar_t e = a;
  detray::detail::set_zero_inverted(e, m);
  for (std::size_t i = 0u; i < scalar_t::size(); ++i) {
    EXPECT_EQ(e[i], i < 3u ? static_cast<value_t>(i + 1u) : 0.f);
  }
  e = a;
  detray::detail::set_zero(e, m);
  for (std::size_t i = 0u; i < scalar_t::size(); ++i) {
    EXPECT_EQ(e[i], i < 3u ? 0.f : static_cast<value_t>(i + 1u));
  }

  // Numeric limits are broadcast
  EXPECT_TRUE(detray::detail::all_of(std::numeric_limits<scalar_t>::epsilon() ==
                                     std::numeric_limits<value_t>::epsilon()));
  EXPECT_TRUE(detray::detail::all_of(std::numeric_limits<scalar_t>::max() ==
                                     std::numeric_limits<value_t>::max()));
  static_assert(std::numeric_limits<scalar_t>::is_specialized);

  // Math functions lane by lane
  const scalar_t v = 0.25f * (idx + 1.f);
  const scalar_t w = -v;
  const scalar_t s_sqrt = math::sqrt(v);
  const scalar_t s_fabs = math::fabs(w);
  const scalar_t s_sin = math::sin(v);
  const scalar_t s_cos = math::cos(v);
  const scalar_t s_atan2 = math::atan2(v, w);
  const scalar_t s_pow = math::pow(v, 2.f);
  const scalar_t s_max = math::max(v, w);
  const scalar_t s_min = math::min(v, w);
  const scalar_t s_cs = math::copysign(v, w);
  const scalar_t s_fma = math::fma(v, v, w);
  const scalar_t s_atanh = math::atanh(0.5f * v);
  const bool_t s_sign = math::signbit(w);
  for (std::size_t i = 0u; i < scalar_t::size(); ++i) {
    EXPECT_NEAR(s_sqrt[i], std::sqrt(v[i]), tol);
    EXPECT_EQ(s_fabs[i], v[i]);
    EXPECT_NEAR(s_sin[i], std::sin(v[i]), tol);
    EXPECT_NEAR(s_cos[i], std::cos(v[i]), tol);
    EXPECT_NEAR(s_atan2[i], std::atan2(v[i], w[i]), tol);
    EXPECT_NEAR(s_pow[i], v[i] * v[i], tol);
    EXPECT_EQ(s_max[i], v[i]);
    EXPECT_EQ(s_min[i], w[i]);
    EXPECT_EQ(s_cs[i], w[i]);
    EXPECT_NEAR(s_fma[i], v[i] * v[i] + w[i], tol);
    EXPECT_NEAR(s_atanh[i], std::atanh(0.5f * v[i]), tol);
    EXPECT_TRUE(s_sign[i]);
  }

  // Single values still resolve to the std functions
  static_assert(std::same_as<decltype(math::sqrt(1.f)), float>);
  static_assert(std::same_as<decltype(math::fabs(-1.)), double>);

  // Approximate comparison
  EXPECT_TRUE(detray::algebra::approx_equal(a, a));
  EXPECT_FALSE(detray::algebra::approx_equal(a, b));

  // Printing
  std::cout << a << std::endl;
  std::cout << m << std::endl;
}

/// This test the vector functions on an SoA (lane bundle) based vector
TEST(detray_algebra_soa, vector) {
  using test_algebra_t = detray::array_soa<value_t>;
  using vector3_v = dvector3D<test_algebra_t>;

  // Value type is a lane bundle of float
  using scalar_t = dscalar<test_algebra_t>;

  static_assert(detray::concepts::scalar<scalar_t>);
  static_assert(detray::concepts::vector<vector3_v>);
  static_assert(detray::concepts::vector3D<vector3_v>);

  // Cast simd scalar to different precisions
  using scalar_f = dscalar<detray::array_soa<float>>;
  using scalar_d = dscalar<detray::array_soa<double>>;
  using scalar_i = dscalar<detray::array_soa<int>>;

  auto s1 = detray::algebra::cast_to<float>(scalar_t(1.f));
  auto s2 = detray::algebra::cast_to<double>(scalar_t(2.f));
  auto s3 = detray::algebra::cast_to<int>(scalar_t(3.f));

  static_assert(std::same_as<decltype(s1), scalar_f>);
  static_assert(std::same_as<decltype(s2), scalar_d>);
  static_assert(std::same_as<decltype(s3), scalar_i>);

  ASSERT_TRUE(detray::detail::all_of(s1 == scalar_f(1.f)));
  ASSERT_TRUE(detray::detail::all_of(s2 == scalar_d(2.)));
  ASSERT_TRUE(detray::detail::all_of(s3 == scalar_i(3)));

  vector3_v a{1.f, 2.f, 3.f};
  vector3_v b{4.f, 5.f, 6.f};

  EXPECT_TRUE(detray::detail::all_of(a[0] == scalar_t(1.f)));
  EXPECT_TRUE(detray::detail::all_of(a[1] == scalar_t(2.f)));
  EXPECT_TRUE(detray::detail::all_of(a[2] == scalar_t(3.f)));

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

  using algebra_f_t = detray::array_soa<float>;
  using algebra_d_t = detray::array_soa<double>;
  using algebra_i_t = detray::array_soa<int>;
  static_assert(std::same_as<decltype(a_cast_f), dvector3D<algebra_f_t>>);
  static_assert(std::same_as<decltype(a_cast_d), dvector3D<algebra_d_t>>);
  static_assert(std::same_as<decltype(a_cast_i), dvector3D<algebra_i_t>>);

  for (int i = 0; i < 3; ++i) {
    EXPECT_TRUE(detray::detail::all_of(a_cast_f[i] ==
                                       detray::algebra::cast_to<float>(a[i])));
    EXPECT_TRUE(detray::detail::all_of(a_cast_d[i] ==
                                       detray::algebra::cast_to<double>(a[i])));
    EXPECT_TRUE(detray::detail::all_of(a_cast_i[i] ==
                                       detray::algebra::cast_to<int>(a[i])));
  }

  // Masked comparison
  auto m = a.compare(a);
  EXPECT_TRUE(detray::detail::all_of(m[0]));
  EXPECT_TRUE(detray::detail::all_of(m[1]));
  EXPECT_TRUE(detray::detail::all_of(m[2]));

  m = a.compare(b);
  EXPECT_FALSE(detray::detail::all_of(m[0]));
  EXPECT_FALSE(detray::detail::all_of(m[1]));
  EXPECT_FALSE(detray::detail::all_of(m[2]));

  // Full comparisons
  EXPECT_TRUE(a == a);
  EXPECT_FALSE(a == b);

  // Addition
  auto v_add = a + b;
  EXPECT_TRUE(detray::detail::all_of(v_add[0] == scalar_t(5.f)));
  EXPECT_TRUE(detray::detail::all_of(v_add[1] == scalar_t(7.f)));
  EXPECT_TRUE(detray::detail::all_of(v_add[2] == scalar_t(9.f)));

  // Subration
  auto v_sub = a - b;
  EXPECT_TRUE(detray::detail::all_of(v_sub[0] == scalar_t(-3.f)));
  EXPECT_TRUE(detray::detail::all_of(v_sub[1] == scalar_t(-3.f)));
  EXPECT_TRUE(detray::detail::all_of(v_sub[2] == scalar_t(-3.f)));

  // Multiplication
  auto v_mul = a * b;
  EXPECT_TRUE(detray::detail::all_of(v_mul[0] == scalar_t(4.f)));
  EXPECT_TRUE(detray::detail::all_of(v_mul[1] == scalar_t(10.f)));
  EXPECT_TRUE(detray::detail::all_of(v_mul[2] == scalar_t(18.f)));

  // Division
  auto v_div = a / b;
  EXPECT_TRUE(detray::detail::all_of(v_div[0] == scalar_t(0.25f)));
  EXPECT_TRUE(detray::detail::all_of(v_div[1] == scalar_t(0.4f)));
  EXPECT_TRUE(detray::detail::all_of(v_div[2] == scalar_t(0.5f)));

  // Scalar multiplication
  auto v_smul = 2.f * b;
  EXPECT_TRUE(detray::detail::all_of(v_smul[0] == scalar_t(8.f)));
  EXPECT_TRUE(detray::detail::all_of(v_smul[1] == scalar_t(10.f)));
  EXPECT_TRUE(detray::detail::all_of(v_smul[2] == scalar_t(12.f)));

  // Lane bundle multiplication
  auto v_bmul = scalar_t(2.f) * b;
  EXPECT_TRUE(detray::detail::all_of(v_bmul[0] == scalar_t(8.f)));
  EXPECT_TRUE(detray::detail::all_of(v_bmul[1] == scalar_t(10.f)));
  EXPECT_TRUE(detray::detail::all_of(v_bmul[2] == scalar_t(12.f)));

  // Expression
  auto v_expr = (b / a) - (2.5f * b) + vector3_v{};
  EXPECT_TRUE(detray::detail::all_of(v_expr[0] == scalar_t(-6.f)));
  EXPECT_TRUE(detray::detail::all_of(v_expr[1] == scalar_t(-10.f)));
  EXPECT_TRUE(detray::detail::all_of(v_expr[2] == scalar_t(-13.f)));

  auto d{vector::dot(a, b)};
  EXPECT_TRUE(detray::detail::all_of(d == scalar_t(32.f)));

  scalar_t norms_a{vector::norm(vector::normalize(a))};
  scalar_t norms_b{vector::norm(vector::normalize(b))};
  for (unsigned int i{0u}; i < norms_a.size(); ++i) {
    EXPECT_NEAR(norms_a[i], 1.f, tol);
    EXPECT_NEAR(norms_b[i], 1.f, tol);
  }

  auto cr{vector::cross(a, b)};
  EXPECT_TRUE(detray::detail::all_of(cr[0] == scalar_t(-3.f)));
  EXPECT_TRUE(detray::detail::all_of(cr[1] == scalar_t(6.f)));
  EXPECT_TRUE(detray::detail::all_of(cr[2] == scalar_t(-3.f)));

  static_assert(std::is_convertible_v<decltype(v_expr), vector3_v>,
                "expression type not convertible");
}

/// This test the getter functions on an SoA (lane bundle) based vector
TEST(detray_algebra_soa, getter) {
  using test_algebra_t = detray::array_soa<value_t>;

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

/// This test an SoA (lane bundle) based affine transform3
TEST(detray_algebra_soa, transform3) {
  using test_algebra_t = detray::array_soa<value_t>;

  // Print the linear algebra types of this backend
  using algebra::operator<<;

  using vector3 = dvector3D<test_algebra_t>;
  using point3 = dpoint3D<test_algebra_t>;
  using scalar_t = dscalar<test_algebra_t>;
  using transform3 = dtransform3D<test_algebra_t>;

  static_assert(detray::concepts::transform3D<transform3>);

  transform3 idty{};

  EXPECT_TRUE(
      detray::detail::all_of(idty(0, 0) == detray::detail::one<scalar_t>()));
  EXPECT_TRUE(
      detray::detail::all_of(idty(1, 0) == detray::detail::zero<scalar_t>()));
  EXPECT_TRUE(
      detray::detail::all_of(idty(2, 0) == detray::detail::zero<scalar_t>()));
  EXPECT_TRUE(
      detray::detail::all_of(idty(0, 1) == detray::detail::zero<scalar_t>()));
  EXPECT_TRUE(
      detray::detail::all_of(idty(1, 1) == detray::detail::one<scalar_t>()));
  EXPECT_TRUE(
      detray::detail::all_of(idty(2, 1) == detray::detail::zero<scalar_t>()));
  EXPECT_TRUE(
      detray::detail::all_of(idty(0, 2) == detray::detail::zero<scalar_t>()));
  EXPECT_TRUE(
      detray::detail::all_of(idty(1, 2) == detray::detail::zero<scalar_t>()));
  EXPECT_TRUE(
      detray::detail::all_of(idty(2, 2) == detray::detail::one<scalar_t>()));
  EXPECT_TRUE(
      detray::detail::all_of(idty(0, 3) == detray::detail::zero<scalar_t>()));
  EXPECT_TRUE(
      detray::detail::all_of(idty(1, 3) == detray::detail::zero<scalar_t>()));
  EXPECT_TRUE(
      detray::detail::all_of(idty(2, 3) == detray::detail::zero<scalar_t>()));

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

  using algebra_f_t = detray::array_soa<float>;
  using algebra_d_t = detray::array_soa<double>;
  using algebra_i_t = detray::array_soa<int>;
  static_assert(std::same_as<decltype(trf1_cast_f), dtransform3D<algebra_f_t>>);
  static_assert(std::same_as<decltype(trf1_cast_d), dtransform3D<algebra_d_t>>);
  static_assert(std::same_as<decltype(trf1_cast_i), dtransform3D<algebra_i_t>>);

  const auto& mat_f = trf1_cast_f.matrix();
  const auto& mat_d = trf1_cast_d.matrix();
  const auto& mat_i = trf1_cast_i.matrix();
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      const auto& elem_ij = trf1.matrix()[i][j];
      EXPECT_TRUE(detray::detail::all_of(
          mat_f[i][j] == detray::algebra::cast_to<float>(elem_ij)));
      EXPECT_TRUE(detray::detail::all_of(
          mat_d[i][j] == detray::algebra::cast_to<double>(elem_ij)));
      EXPECT_TRUE(detray::detail::all_of(
          mat_i[i][j] == detray::algebra::cast_to<int>(elem_ij)));
    }
  }

  EXPECT_TRUE(detray::detail::all_of(trf2(0, 0) == x[0]));
  EXPECT_TRUE(detray::detail::all_of(trf2(1, 0) == x[1]));
  EXPECT_TRUE(detray::detail::all_of(trf2(2, 0) == x[2]));
  EXPECT_TRUE(detray::detail::all_of(trf2(0, 1) == y[0]));
  EXPECT_TRUE(detray::detail::all_of(trf2(1, 1) == y[1]));
  EXPECT_TRUE(detray::detail::all_of(trf2(2, 1) == y[2]));
  EXPECT_TRUE(detray::detail::all_of(trf2(0, 2) == z[0]));
  EXPECT_TRUE(detray::detail::all_of(trf2(1, 2) == z[1]));
  EXPECT_TRUE(detray::detail::all_of(trf2(2, 2) == z[2]));
  EXPECT_TRUE(detray::detail::all_of(trf2(0, 3) ==
                                     2.f * detray::detail::one<scalar_t>()));
  EXPECT_TRUE(detray::detail::all_of(trf2(1, 3) ==
                                     3.f * detray::detail::one<scalar_t>()));
  EXPECT_TRUE(detray::detail::all_of(trf2(2, 3) ==
                                     4.f * detray::detail::one<scalar_t>()));

  // Check that local origin translates into global translation
  point3 lzero = {0.f, 0.f, 0.f};
  point3 gzero = trf2.point_to_global(lzero);
  EXPECT_TRUE(detray::detail::all_of(gzero[0] == t[0]));
  EXPECT_TRUE(detray::detail::all_of(gzero[1] == t[1]));
  EXPECT_TRUE(detray::detail::all_of(gzero[2] == t[2]));

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

/// This test an SoA (lane bundle) based 2x3 matrix
TEST(detray_algebra_soa, matrix3) {
  using test_algebra_t = detray::array_soa<value_t>;
  using matrix_2x3_t = dmatrix<test_algebra_t, 2, 3>;

  // Test type traits
  static_assert(
      std::is_same_v<detray::traits::index_t<matrix_2x3_t>, std::size_t>);
  static_assert(std::is_same_v<detray::traits::value_t<matrix_2x3_t>, value_t>);
  static_assert(std::is_same_v<detray::traits::scalar_t<matrix_2x3_t>,
                               dscalar<test_algebra_t>>);
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

/// This test an SoA (lane bundle) based 6x4 matrix
TEST(detray_algebra_soa, matrix64) {
  using test_algebra_t = detray::array_soa<value_t>;

  // Create the matrix.
  using matrix_6x4_t = dmatrix<test_algebra_t, 6, 4>;
  using scalar_t = dscalar<test_algebra_t>;

  matrix_6x4_t m;

  // Test type traits
  static_assert(
      std::is_same_v<detray::traits::index_t<matrix_6x4_t>, std::size_t>);
  static_assert(std::is_same_v<detray::traits::value_t<matrix_6x4_t>, value_t>);
  static_assert(
      std::is_same_v<detray::traits::scalar_t<matrix_6x4_t>, scalar_t>);

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

  using algebra_f_t = detray::array_soa<float>;
  using algebra_d_t = detray::array_soa<double>;
  using algebra_i_t = detray::array_soa<int>;
  static_assert(std::same_as<decltype(m_cast_f), dmatrix<algebra_f_t, 6, 4>>);
  static_assert(std::same_as<decltype(m_cast_d), dmatrix<algebra_d_t, 6, 4>>);
  static_assert(std::same_as<decltype(m_cast_i), dmatrix<algebra_i_t, 6, 4>>);

  // Matrix products and transpose on lane bundles
  using matrix_4x6_t = dmatrix<test_algebra_t, 4, 6>;
  const matrix_4x6_t I64_T = detray::matrix::transpose(I64);
  for (std::size_t i = 0u; i < 4u; ++i) {
    for (std::size_t j = 0u; j < 6u; ++j) {
      EXPECT_TRUE(detray::detail::all_of(getter::element(I64_T, i, j) ==
                                         getter::element(I64, j, i)));
    }
  }
  // Only the first four diagonal entries of the product are one
  const auto I66 = I64 * I64_T;
  for (std::size_t i = 0u; i < 6u; ++i) {
    for (std::size_t j = 0u; j < 6u; ++j) {
      const value_t ref = (i == j && i < 4u) ? 1.f : 0.f;
      EXPECT_TRUE(detray::detail::all_of(getter::element(I66, i, j) == ref))
          << i << j;
    }
  }
}

/// This tests the SoA branch of the quadratic equation solver
TEST(detray_algebra_soa, quadratic_equation) {
  using test_algebra_t = detray::array_soa<value_t>;
  using scalar_t = dscalar<test_algebra_t>;

  // Lane 2: one solution (double root), lane 3: two solutions with a
  // negative leading coefficient, the rest: two solutions
  scalar_t a{1.f};
  a[3] = -1.f;
  scalar_t b{-3.f};
  b[2] = 2.f;
  scalar_t c{2.f};
  c[2] = 1.f;
  c[3] = -2.f;

  const detail::quadratic_equation<scalar_t> qe{a, b, c};

  for (std::size_t i = 0u; i < scalar_t::size(); ++i) {
    const detail::quadratic_equation<value_t> qe_ref{a[i], b[i], c[i]};

    EXPECT_EQ(static_cast<int>(qe.solutions()[i]), qe_ref.solutions()) << i;
    EXPECT_NEAR(qe.smaller()[i], qe_ref.smaller(), tol) << i;
    if (qe_ref.solutions() == 2) {
      EXPECT_NEAR(qe.larger()[i], qe_ref.larger(), tol) << i;
    }
  }

  // All lanes linear: early exit
  const detail::quadratic_equation<scalar_t> qe_lin{
      scalar_t{0.f}, scalar_t{2.f}, scalar_t{-4.f}};
  EXPECT_TRUE(detray::detail::all_of(qe_lin.solutions() == 1.f));
  EXPECT_TRUE(detray::detail::all_of(qe_lin.smaller() == 2.f));
}

/// This tests that the SoA ray intersectors compile and work with the plugin
TEST(detray_algebra_soa, ray_plane_intersection) {
  using test_algebra_t = detray::array_soa<value_t>;
  using scalar_t = dscalar<test_algebra_t>;
  using point3_t = dpoint3D<test_algebra_t>;
  using vector3_t = dvector3D<test_algebra_t>;
  using transform3_t = dtransform3D<test_algebra_t>;

  enum class mask_id : unsigned int { e_rectangle2D = 0 };
  enum class material_id : unsigned int { e_slab = 0 };
  using mask_link_t = dtyped_index<mask_id, dindex>;
  using material_link_t = dtyped_index<material_id, dindex>;
  using surface_desc_t = surface_descriptor<mask_link_t, material_link_t>;
  using mask_t = mask<rectangle2D, test_algebra_t, std::uint_least16_t>;

  const surface_desc_t sf(0u, mask_link_t{mask_id::e_rectangle2D, 0u},
                          material_link_t{material_id::e_slab, 0u}, 0u,
                          surface_id::e_sensitive);
  // Half lengths 5 x 5: the lanes with |x| > 5 are outside
  const mask_t rect{0u, 5.f, 5.f};

  // One plane per lane, shifted in x by 3 * lane and in z by 10 + lane
  const scalar_t idx = detray::detail::iota<scalar_t>();
  const point3_t t{3.f * idx, 0.f, 10.f + idx};
  const vector3_t z{0.f, 0.f, 1.f};
  const vector3_t x{1.f, 0.f, 0.f};
  const transform3_t trf(t, z, x);

  // Ray along z from the origin (broadcast to all lanes)
  const detail::ray<test_algebra_t> r{point3_t{0.f, 0.f, 0.f}, 0.f,
                                      vector3_t{0.f, 0.f, 1.f}, 0.f};

  const auto is =
      ray_intersector<rectangle2D, test_algebra_t>{}(r, sf, rect, trf);

  EXPECT_TRUE(is.is_inside());
  EXPECT_FALSE(is.is_outside());
  EXPECT_TRUE(detail::all_of(is.is_along()));

  for (std::size_t i = 0u; i < scalar_t::size(); ++i) {
    EXPECT_NEAR(is.path()[i], 10.f + static_cast<value_t>(i), tol) << i;

    const auto status = static_cast<intersection::status>(is.status()[i]);
    if (3.f * static_cast<value_t>(i) <= 5.f) {
      EXPECT_EQ(status, intersection::status::e_inside) << i;
    } else {
      EXPECT_EQ(status, intersection::status::e_outside) << i;
    }
  }
}
