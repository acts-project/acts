// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Test include(s).
#include "detray/test/cpu/algebra_fixture.hpp"

// Project include(s).
#include "detray/algebra/utils/approximately_equal.hpp"
#include "detray/algebra/utils/casts.hpp"
#include "detray/algebra/utils/print.hpp"

// GoogleTest include(s).
#include <gtest/gtest.h>

namespace detray::test {

// This defines the local frame test suite
TEST_F(detray_algebra, vector2D) {
  // Construction
  point2 vA{0.f, 1.f};
  ASSERT_EQ(vA[0], 0.f);
  ASSERT_EQ(vA[1], 1.f);

  // Test printing
  std::cout << vA << std::endl;

  // Test comparison
  ASSERT_TRUE(detray::algebra::approx_equal(vA, vA));
  ASSERT_TRUE(detray::algebra::approx_equal(vA, vA, this->epsilon()));

  scalar rel_err{1.f + 10.f * this->epsilon()};
  point2 vA_err = rel_err * vA;
  ASSERT_TRUE(
      detray::algebra::approx_equal(vA, vA_err, 11.f * this->epsilon()));
  ASSERT_FALSE(
      detray::algebra::approx_equal(vA, vA_err, 9.f * this->epsilon()));

  rel_err = 1.f + 17.f * this->epsilon();
  vA_err = rel_err * vA;
  ASSERT_TRUE(
      detray::algebra::approx_equal(vA, vA_err, 18.f * this->epsilon()));
  ASSERT_FALSE(
      detray::algebra::approx_equal(vA, vA_err, 16.f * this->epsilon()));
  // Cast to (different) precision
  const auto vA_cast_f = detray::algebra::cast_to<float>(vA);

  for (index i = 0; i < 2; ++i) {
    auto elem_i = vA_cast_f[i];

    static_assert(std::same_as<decltype(elem_i), float>);
    ASSERT_FLOAT_EQ(elem_i, static_cast<float>(vA[i]));
  }

  const auto vA_cast_d = detray::algebra::cast_to<double>(vA);

  for (index i = 0; i < 2; ++i) {
    auto elem_i = vA_cast_d[i];

    static_assert(std::same_as<decltype(elem_i), double>);
    ASSERT_DOUBLE_EQ(elem_i, static_cast<double>(vA[i]));
  }

  const auto vA_cast_i = detray::algebra::cast_to<int>(vA);

  for (index i = 0; i < 2; ++i) {
    auto elem_i = vA_cast_i[i];

    static_assert(std::same_as<decltype(elem_i), int>);
    ASSERT_EQ(elem_i, static_cast<int>(vA[i]));
  }

  // Assignment
  point2 vB = vA;
  ASSERT_EQ(vB[0], 0.f);
  ASSERT_EQ(vB[1], 1.f);

  // Addition
  point2 vC = vA + vB;
  ASSERT_EQ(vC[0], 0.f);
  ASSERT_EQ(vC[1], 2.f);

  // Multiplication by scalar
  point2 vC2 = vC * 2.f;
  ASSERT_EQ(vC2[0], 0.f);
  ASSERT_EQ(vC2[1], 4.f);

  // Cast operations to phi, theta, eta, perp
  vector2 vD{1.f, 1.f};
  scalar phi = detray::vector::phi(vD);
  ASSERT_NEAR(phi, constant<scalar>::pi_4, this->epsilon());

  scalar perp = detray::vector::perp(vD);
  ASSERT_NEAR(perp, std::sqrt(2.), this->epsilon());

  scalar norm = detray::vector::norm(vD);
  ASSERT_NEAR(norm, std::sqrt(2.), this->epsilon());

  vector2 vDnorm = detray::vector::normalize(vD);
  ASSERT_NEAR(vDnorm[0], 1. / std::sqrt(2.), this->epsilon());
  ASSERT_NEAR(vDnorm[1], 1. / std::sqrt(2.), this->epsilon());

  // Test comparison
  point2 vE{100.f, 462809.f};

  ASSERT_TRUE(detray::algebra::approx_equal(vE, vE));
  ASSERT_TRUE(detray::algebra::approx_equal(vE, vE, this->epsilon()));

  rel_err = 1.f + 10.f * this->epsilon();
  point2 vE_err = rel_err * vE;
  ASSERT_TRUE(
      detray::algebra::approx_equal(vE, vE_err, 11.f * this->epsilon()));
  ASSERT_FALSE(
      detray::algebra::approx_equal(vE, vE_err, 9.f * this->epsilon()));

  point2 vE_abs_err{100.00001f, 462809.05f};
  ASSERT_TRUE(detray::algebra::approx_equal(vE, vE_abs_err, 0.00001f));
  ASSERT_FALSE(detray::algebra::approx_equal(vE, vE_abs_err, this->epsilon()));

  rel_err = 1.f + 17.f * this->epsilon();
  vE_err = rel_err * vE;
  ASSERT_TRUE(
      detray::algebra::approx_equal(vE, vE_err, 18.f * this->epsilon()));
  ASSERT_FALSE(
      detray::algebra::approx_equal(vE, vE_err, 16.f * this->epsilon()));
}

// This defines the vector3 test suite
TEST_F(detray_algebra, vector3D) {
  // Test concepts
  static_assert(detray::concepts::scalar<scalar>);
  static_assert(!detray::concepts::vector<scalar>);

  static_assert(!detray::concepts::scalar<vector3>);
  static_assert(detray::concepts::vector<vector3>);
  static_assert(detray::concepts::vector3D<vector3>);
  static_assert(!detray::concepts::vector2D<vector3>);

  static_assert(!detray::concepts::scalar<vector2>);
  static_assert(detray::concepts::vector<vector2>);
  static_assert(detray::concepts::vector2D<vector2>);
  static_assert(!detray::concepts::vector3D<vector2>);

  // Construction
  vector3 vA{0.f, 1.f, 2.f};
  ASSERT_EQ(vA[0], 0.f);
  ASSERT_EQ(vA[1], 1.f);
  ASSERT_EQ(vA[2], 2.f);

  // Test printing
  std::cout << vA << std::endl;

  // Cast to (different) precision
  const auto vA_cast_f = detray::algebra::cast_to<float>(vA);

  for (index i = 0; i < 3; ++i) {
    auto elem_i = vA_cast_f[i];

    static_assert(std::same_as<decltype(elem_i), float>);
    ASSERT_FLOAT_EQ(elem_i, static_cast<float>(vA[i]));
  }

  const auto vA_cast_d = detray::algebra::cast_to<double>(vA);

  for (index i = 0; i < 3; ++i) {
    auto elem_i = vA_cast_d[i];

    static_assert(std::same_as<decltype(elem_i), double>);
    ASSERT_DOUBLE_EQ(elem_i, static_cast<double>(vA[i]));
  }

  const auto vA_cast_i = detray::algebra::cast_to<int>(vA);

  for (index i = 0; i < 3; ++i) {
    auto elem_i = vA_cast_i[i];

    static_assert(std::same_as<decltype(elem_i), int>);
    ASSERT_EQ(elem_i, static_cast<int>(vA[i]));
  }

  // Assignment
  vector3 vB = vA;
  ASSERT_EQ(vB[0], 0.f);
  ASSERT_EQ(vB[1], 1.f);
  ASSERT_EQ(vB[2], 2.f);

  // Addition
  vector3 vC = vA + vB;
  ASSERT_EQ(vC[0], 0.f);
  ASSERT_EQ(vC[1], 2.f);
  ASSERT_EQ(vC[2], 4.f);

  // Multiplication by scalar
  vector3 vC2 = vC * 2.0f;
  ASSERT_EQ(vC2[0], 0.f);
  ASSERT_EQ(vC2[1], 4.f);
  ASSERT_EQ(vC2[2], 8.f);

  // Cast operations to phi, theta, eta, perp
  vector3 vD{1.f, 1.f, 1.f};
  scalar phi = detray::vector::phi(vD);
  ASSERT_NEAR(phi, constant<scalar>::pi_4, this->epsilon());

  scalar theta = detray::vector::theta(vD);
  ASSERT_NEAR(theta, std::atan2(std::sqrt(2.), 1.), this->epsilon());

  scalar eta = detray::vector::eta(vD);
  ASSERT_NEAR(eta, 0.65847891569137573, this->tolerance());

  scalar perp = detray::vector::perp(vD);
  ASSERT_NEAR(perp, std::sqrt(2.), this->epsilon());

  scalar norm = detray::vector::norm(vD);
  ASSERT_NEAR(norm, std::sqrt(3.), this->epsilon());
}

// This defines the vector operation test suite
TEST_F(detray_algebra, element_getter) {
  vector3 v3{1.f, 1.f, 1.f};

  // Normalization
  vector3 v3n = detray::vector::normalize(v3);
  ASSERT_NEAR(v3n[0], 1. / std::sqrt(3.), this->epsilon());
  ASSERT_NEAR(v3n[1], 1. / std::sqrt(3.), this->epsilon());
  ASSERT_NEAR(v3n[2], 1. / std::sqrt(3.), this->epsilon());

  // Cross product
  vector3 z = detray::vector::normalize(vector3{3.f, 2.f, 1.f});
  vector3 x = detray::vector::normalize(vector3{2.f, -3.f, 0.f});
  vector3 y = detray::vector::cross(z, detray::vector::normalize(x));

  // Check with dot product
  ASSERT_NEAR(detray::vector::dot(x, y), 0.f, this->epsilon());
  ASSERT_NEAR(detray::vector::dot(y, z), 0.f, this->epsilon());
  ASSERT_NEAR(detray::vector::dot(z, x), 0.f, this->epsilon());
}

// This test checks to see if the `dot` function can handle when one of its
// operands is a sum or difference of two vectors. We also test to see if the
// `dot` function plays nicely with scalar multiplication (just to be safe).
TEST_F(detray_algebra, dot_product_with_ops) {
  vector3 v1{1.f, 2.f, 3.f};
  vector3 v2{3.f, 4.f, 5.f};

  ASSERT_NEAR(detray::vector::dot(v1 + v2, v2), 76.f, this->epsilon());
  ASSERT_NEAR(detray::vector::dot(v1, v2 - v1), 12.f, this->epsilon());
  ASSERT_NEAR(detray::vector::dot(v1 + v2, v1 - v2), -36.f, this->epsilon());
  ASSERT_NEAR(detray::vector::dot(v1 + v2, 2 * v2), 152.f, this->epsilon());
}

// This test checks to see if the `cross` function can handle when one of its
// operands is a sum or difference of two vectors. We also test to see if the
// `cross` function plays nicely with scalar multiplication (just to be safe).
TEST_F(detray_algebra, cross_product_add_sub) {
  vector3 v1{1.f, 2.f, 3.f};
  vector3 v2{3.f, 4.f, 5.f};
  vector3 v3{-6.f, 7.f, -9.f};

  vector3 v = detray::vector::cross(v1 + v2, v3);
  vector3 ans{-110.f, -12.f, 64.f};

  ASSERT_NEAR(v[0], ans[0], this->epsilon());
  ASSERT_NEAR(v[1], ans[1], this->epsilon());
  ASSERT_NEAR(v[2], ans[2], this->epsilon());

  v = detray::vector::cross(v3 - 2 * v1, 3 * (v1 + v2));
  ans = {342.f, 12.f, -180.f};

  ASSERT_NEAR(v[0], ans[0], this->epsilon());
  ASSERT_NEAR(v[1], ans[1], this->epsilon());
  ASSERT_NEAR(v[2], ans[2], this->epsilon());
}

}  // namespace detray::test
