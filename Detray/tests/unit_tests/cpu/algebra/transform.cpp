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

// Test include(s).
#include "detray/test/cpu/algebra_fixture.hpp"

// Project include(s).
#include "detray/algebra/utils/approximately_equal.hpp"
#include "detray/algebra/utils/casts.hpp"
#include "detray/algebra/utils/print.hpp"

// GoogleTest include(s).
#include <gtest/gtest.h>

// System include(s)
#include <array>
#include <iostream>

namespace detray::test {

// This defines the transform3 test suite
TEST_F(detray_algebra, transform3D) {
  static_assert(detray::concepts::transform3D<transform3>);

  // Preparation work
  vector3 z = detray::vector::normalize(vector3{3.f, 2.f, 1.f});
  vector3 x = detray::vector::normalize(vector3{2.f, -3.f, 0.f});
  vector3 y = detray::vector::cross(z, x);
  point3 t = {2.f, 3.f, 4.f};

  // Test constructor from t, z, x
  transform3 trf1(t, z, x);
  ASSERT_TRUE(trf1 == trf1);
  transform3 trf2;
  trf2 = trf1;

  // Test printing
  std::cout << trf1 << std::endl;

  // Test comparison
  ASSERT_TRUE(detray::algebra::approx_equal(trf1, trf1));
  ASSERT_TRUE(detray::algebra::approx_equal(trf1, trf1, this->epsilon()));

  scalar rel_err{1.f + 10.f * this->epsilon()};
  transform3 trf1_err(rel_err * t, rel_err * z, rel_err * x);
  ASSERT_TRUE(
      detray::algebra::approx_equal(trf1, trf1_err, 200.f * this->epsilon()));
  ASSERT_FALSE(
      detray::algebra::approx_equal(trf1, trf1_err, 10.f * this->epsilon()));
  // Cast to (different) precision
  const auto trf1_cast_f = detray::algebra::cast_to<float>(trf1);
  const auto& mat_f = trf1_cast_f.matrix();
  const auto& mat_inv_f = trf1_cast_f.matrix_inverse();

  for (std::size_t j = 0u; j < 3u; ++j) {
    for (std::size_t i = 0u; i < 2u; ++i) {
      auto elem_i = detray::getter::element(mat_f, i, j);
      auto elem_inv_i = detray::getter::element(mat_inv_f, i, j);

      static_assert(std::same_as<decltype(elem_i), float>);
      static_assert(std::same_as<decltype(elem_inv_i), float>);
      ASSERT_FLOAT_EQ(elem_i,
                      static_cast<float>(detray::getter::element(mat_f, i, j)));
      ASSERT_FLOAT_EQ(elem_inv_i, static_cast<float>(detray::getter::element(
                                      mat_inv_f, i, j)));
    }
  }

  const auto trf1_cast_d = detray::algebra::cast_to<double>(trf1);
  const auto& mat_d = trf1_cast_d.matrix();
  const auto& mat_inv_d = trf1_cast_d.matrix_inverse();

  for (std::size_t j = 0u; j < 3u; ++j) {
    for (std::size_t i = 0u; i < 2u; ++i) {
      auto elem_i = detray::getter::element(mat_d, i, j);
      auto elem_inv_i = detray::getter::element(mat_inv_d, i, j);

      static_assert(std::same_as<decltype(elem_i), double>);
      static_assert(std::same_as<decltype(elem_inv_i), double>);
      ASSERT_DOUBLE_EQ(
          elem_i, static_cast<double>(detray::getter::element(mat_d, i, j)));
      ASSERT_DOUBLE_EQ(elem_inv_i, static_cast<double>(detray::getter::element(
                                       mat_inv_d, i, j)));
    }
  }

  const auto trf1_cast_i = detray::algebra::cast_to<int>(trf1);
  const auto& mat_i = trf1_cast_i.matrix();
  const auto& mat_inv_i = trf1_cast_i.matrix_inverse();

  for (std::size_t j = 0u; j < 3u; ++j) {
    for (std::size_t i = 0u; i < 2u; ++i) {
      auto elem_i = detray::getter::element(mat_i, i, j);
      auto elem_inv_i = detray::getter::element(mat_inv_i, i, j);

      static_assert(std::same_as<decltype(elem_i), int>);
      static_assert(std::same_as<decltype(elem_inv_i), int>);
      ASSERT_EQ(elem_i, static_cast<int>(detray::getter::element(mat_i, i, j)));
      ASSERT_EQ(elem_inv_i,
                static_cast<int>(detray::getter::element(mat_inv_i, i, j)));
    }
  }

  const auto rot = trf2.rotation();
  ASSERT_NEAR(detray::getter::element(rot, 0, 0), x[0], this->epsilon());
  ASSERT_NEAR(detray::getter::element(rot, 1, 0), x[1], this->epsilon());
  ASSERT_NEAR(detray::getter::element(rot, 2, 0), x[2], this->epsilon());
  ASSERT_NEAR(detray::getter::element(rot, 0, 1), y[0], this->epsilon());
  ASSERT_NEAR(detray::getter::element(rot, 1, 1), y[1], this->epsilon());
  ASSERT_NEAR(detray::getter::element(rot, 2, 1), y[2], this->epsilon());
  ASSERT_NEAR(detray::getter::element(rot, 0, 2), z[0], this->epsilon());
  ASSERT_NEAR(detray::getter::element(rot, 1, 2), z[1], this->epsilon());
  ASSERT_NEAR(detray::getter::element(rot, 2, 2), z[2], this->epsilon());

  auto trn = trf2.translation();
  ASSERT_NEAR(trn[0], 2.f, this->epsilon());
  ASSERT_NEAR(trn[1], 3.f, this->epsilon());
  ASSERT_NEAR(trn[2], 4.f, this->epsilon());

  // Test constructor from matrix
  auto m44 = trf2.matrix();
  transform3 trfm(m44);

  // Make sure that detray::getter:vector can be called.
  [[maybe_unused]] vector3
      test_vector;  // we need to declare a variable in order to use the
                    // [[maybe_unused]] attribute here

  test_vector = detray::getter::vector<3>(m44, 0, 0);
  test_vector = detray::getter::vector<3>(m44, 0, 1);
  test_vector = detray::getter::vector<3>(m44, 0, 2);

  // Test constructor from inverse matrix
  auto m44_inv = trf2.matrix_inverse();

  // Make sure that detray::getter:vector can be called.
  test_vector = detray::getter::vector<3>(m44_inv, 0, 0);
  test_vector = detray::getter::vector<3>(m44_inv, 0, 1);
  test_vector = detray::getter::vector<3>(m44_inv, 0, 2);

  // Re-evaluate rot and trn
  auto rotm = trfm.rotation();
  ASSERT_NEAR(detray::getter::element(rotm, 0, 0), x[0], this->epsilon());
  ASSERT_NEAR(detray::getter::element(rotm, 1, 0), x[1], this->epsilon());
  ASSERT_NEAR(detray::getter::element(rotm, 2, 0), x[2], this->epsilon());
  ASSERT_NEAR(detray::getter::element(rotm, 0, 1), y[0], this->epsilon());
  ASSERT_NEAR(detray::getter::element(rotm, 1, 1), y[1], this->epsilon());
  ASSERT_NEAR(detray::getter::element(rotm, 2, 1), y[2], this->epsilon());
  ASSERT_NEAR(detray::getter::element(rotm, 0, 2), z[0], this->epsilon());
  ASSERT_NEAR(detray::getter::element(rotm, 1, 2), z[1], this->epsilon());
  ASSERT_NEAR(detray::getter::element(rotm, 2, 2), z[2], this->epsilon());

  auto trnm = trfm.translation();
  ASSERT_NEAR(trnm[0], 2.f, this->epsilon());
  ASSERT_NEAR(trnm[1], 3.f, this->epsilon());
  ASSERT_NEAR(trnm[2], 4.f, this->epsilon());

  // Check a construction from an array[16]
  std::array<scalar, 16> matray_helper = {1.f, 0.f, 0.f, 0.f, 0.f, 1.f,
                                          0.f, 0.f, 0.f, 0.f, 1.f, 0.f,
                                          0.f, 0.f, 0.f, 0.f};
  typename transform3::template array_type<16> matray;
  for (unsigned int i = 0u; i < 16u; ++i) {
    matray[i] = matray_helper[i];
  }
  transform3 trfma(matray);

  // Re-evaluate rot and trn
  auto rotma = trfma.rotation();
  ASSERT_NEAR(detray::getter::element(rotma, 0, 0), 1.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(rotma, 1, 0), 0.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(rotma, 2, 0), 0.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(rotma, 0, 1), 0.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(rotma, 1, 1), 1.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(rotma, 2, 1), 0.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(rotma, 0, 2), 0.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(rotma, 1, 2), 0.f, this->epsilon());
  ASSERT_NEAR(detray::getter::element(rotma, 2, 2), 1.f, this->epsilon());

  auto trnma = trfma.translation();
  ASSERT_NEAR(trnma[0], 0.f, this->epsilon());
  ASSERT_NEAR(trnma[1], 0.f, this->epsilon());
  ASSERT_NEAR(trnma[2], 0.f, this->epsilon());
}

// This test global coordinate transforms
TEST_F(detray_algebra, global_transformations) {
  // Preparation work
  vector3 z = detray::vector::normalize(vector3{3.f, 2.f, 1.f});
  vector3 x = detray::vector::normalize(vector3{2.f, -3.f, 0.f});
  [[maybe_unused]] vector3 y = detray::vector::cross(z, x);
  point3 t = {2.f, 3.f, 4.f};
  transform3 trf(t, z, x);

  // Check that local origin translates into global translation
  point3 lzero = {0.f, 0.f, 0.f};
  point3 gzero = trf.point_to_global(lzero);
  ASSERT_NEAR(gzero[0], t[0], this->epsilon());
  ASSERT_NEAR(gzero[1], t[1], this->epsilon());
  ASSERT_NEAR(gzero[2], t[2], this->epsilon());

  // Check a round trip for point
  point3 lpoint = {3.f, 4.f, 5.f};
  point3 gpoint = trf.point_to_global(lpoint);
  point3 lpoint_r = trf.point_to_local(gpoint);
  ASSERT_NEAR(lpoint[0], lpoint_r[0], this->tolerance());
  ASSERT_NEAR(lpoint[1], lpoint_r[1], this->tolerance());
  ASSERT_NEAR(lpoint[2], lpoint_r[2], this->tolerance());

  // Check a point versus vector transform
  // vector should not change if transformed by a pure translation
  transform3 ttrf(t);

  vector3 gvector = {1.f, 1.f, 1.f};
  vector3 lvector = ttrf.vector_to_local(gvector);
  ASSERT_NEAR(gvector[0], lvector[0], this->tolerance());
  ASSERT_NEAR(gvector[1], lvector[1], this->tolerance());
  ASSERT_NEAR(gvector[2], lvector[2], this->tolerance());

  // Check a round trip for vector
  vector3 lvectorB = {7.f, 8.f, 9.f};
  vector3 gvectorB = trf.vector_to_local(lvectorB);
  vector3 lvectorC = trf.vector_to_global(gvectorB);
  ASSERT_NEAR(lvectorB[0], lvectorC[0], this->tolerance());
  ASSERT_NEAR(lvectorB[1], lvectorC[1], this->tolerance());
  ASSERT_NEAR(lvectorB[2], lvectorC[2], this->tolerance());
}

}  // namespace detray::test
