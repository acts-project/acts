// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Detray test include(s)
#include "detray/test/cpu/algebra_fixture.hpp"

// GTest include(s)
#include <gtest/gtest.h>

namespace detray::test {

TEST_F(detray_algebra, column_wise_cross) {
  using matrix3_type = matrix<3, 3>;

  auto P = detray::matrix::zero<matrix3_type>();

  getter::element(P, 0u, 0u) = 0.f;
  getter::element(P, 0u, 1u) = 1.f;
  getter::element(P, 0u, 2u) = 2.f;
  getter::element(P, 1u, 0u) = 3.f;
  getter::element(P, 1u, 1u) = 4.f;
  getter::element(P, 1u, 2u) = 5.f;
  getter::element(P, 2u, 0u) = 6.f;
  getter::element(P, 2u, 1u) = 7.f;
  getter::element(P, 2u, 2u) = 8.f;

  const vector3 u{1.f, 2.f, 3.f};

  const matrix3_type Q = detray::matrix::column_wise_cross(P, u);

  EXPECT_NEAR(getter::element(Q, 0u, 0u), -3.f, this->tolerance());
  EXPECT_NEAR(getter::element(Q, 1u, 0u), 6.f, this->tolerance());
  EXPECT_NEAR(getter::element(Q, 2u, 0u), -3.f, this->tolerance());
  EXPECT_NEAR(getter::element(Q, 0u, 1u), -2.f, this->tolerance());
  EXPECT_NEAR(getter::element(Q, 1u, 1u), 4.f, this->tolerance());
  EXPECT_NEAR(getter::element(Q, 2u, 1u), -2.f, this->tolerance());
  EXPECT_NEAR(getter::element(Q, 0u, 2u), -1.f, this->tolerance());
  EXPECT_NEAR(getter::element(Q, 1u, 2u), 2.f, this->tolerance());
  EXPECT_NEAR(getter::element(Q, 2u, 2u), -1.f, this->tolerance());
}

TEST_F(detray_algebra, column_wise_multiply) {
  using matrix3_type = matrix<3, 3>;

  auto P = detray::matrix::zero<matrix3_type>();

  getter::element(P, 0u, 0u) = 0.f;
  getter::element(P, 0u, 1u) = 1.f;
  getter::element(P, 0u, 2u) = 2.f;
  getter::element(P, 1u, 0u) = 3.f;
  getter::element(P, 1u, 1u) = 4.f;
  getter::element(P, 1u, 2u) = 5.f;
  getter::element(P, 2u, 0u) = 6.f;
  getter::element(P, 2u, 1u) = 7.f;
  getter::element(P, 2u, 2u) = 8.f;

  const vector3 u{1.f, 2.f, 3.f};

  const matrix3_type Q = detray::matrix::column_wise_multiply(P, u);

  EXPECT_NEAR(getter::element(Q, 0u, 0u), 0.f, this->tolerance());
  EXPECT_NEAR(getter::element(Q, 0u, 1u), 1.f, this->tolerance());
  EXPECT_NEAR(getter::element(Q, 0u, 2u), 2.f, this->tolerance());
  EXPECT_NEAR(getter::element(Q, 1u, 0u), 6.f, this->tolerance());
  EXPECT_NEAR(getter::element(Q, 1u, 1u), 8.f, this->tolerance());
  EXPECT_NEAR(getter::element(Q, 1u, 2u), 10.f, this->tolerance());
  EXPECT_NEAR(getter::element(Q, 2u, 0u), 18.f, this->tolerance());
  EXPECT_NEAR(getter::element(Q, 2u, 1u), 21.f, this->tolerance());
  EXPECT_NEAR(getter::element(Q, 2u, 2u), 24.f, this->tolerance());
}

TEST_F(detray_algebra, cross_matrix) {
  using matrix3_type = matrix<3, 3>;

  const vector3 u{1.f, 2.f, 3.f};
  const vector3 v{3.f, 4.f, 5.f};

  const matrix3_type u_cross = detray::matrix::cross_matrix(u);
  const matrix3_type v_cross = detray::matrix::cross_matrix(v);

  EXPECT_NEAR(getter::element(u_cross, 0u, 0u), 0.f, this->tolerance());
  EXPECT_NEAR(getter::element(u_cross, 0u, 1u), -3.f, this->tolerance());
  EXPECT_NEAR(getter::element(u_cross, 0u, 2u), 2.f, this->tolerance());
  EXPECT_NEAR(getter::element(u_cross, 1u, 0u), 3.f, this->tolerance());
  EXPECT_NEAR(getter::element(u_cross, 1u, 1u), 0.f, this->tolerance());
  EXPECT_NEAR(getter::element(u_cross, 1u, 2u), -1.f, this->tolerance());
  EXPECT_NEAR(getter::element(u_cross, 2u, 0u), -2.f, this->tolerance());
  EXPECT_NEAR(getter::element(u_cross, 2u, 1u), 1.f, this->tolerance());
  EXPECT_NEAR(getter::element(u_cross, 2u, 2u), 0.f, this->tolerance());

  // [u]_cross * v = [v]_cross^T * u
  const vector3 u_cross_v = u_cross * v;
  const vector3 v_cross_u = detray::matrix::transpose(v_cross) * u;

  EXPECT_NEAR(getter::element(u_cross_v, 0u), -2.f, this->tolerance());
  EXPECT_NEAR(getter::element(u_cross_v, 1u), 4.f, this->tolerance());
  EXPECT_NEAR(getter::element(u_cross_v, 2u), -2.f, this->tolerance());
  EXPECT_NEAR(getter::element(u_cross_v, 0u), getter::element(v_cross_u, 0u),
              this->tolerance());
  EXPECT_NEAR(getter::element(u_cross_v, 1u), getter::element(v_cross_u, 1u),
              this->tolerance());
  EXPECT_NEAR(getter::element(u_cross_v, 2u), getter::element(v_cross_u, 2u),
              this->tolerance());
}

TEST_F(detray_algebra, outer_product) {
  using matrix3_type = matrix<3, 3>;

  const vector3 u{1.f, 2.f, 3.f};
  const vector3 v{3.f, 4.f, 5.f};

  const matrix3_type m33 = detray::matrix::outer_product(u, v);

  EXPECT_NEAR(getter::element(m33, 0u, 0u), 3.f, this->tolerance());
  EXPECT_NEAR(getter::element(m33, 0u, 1u), 4.f, this->tolerance());
  EXPECT_NEAR(getter::element(m33, 0u, 2u), 5.f, this->tolerance());
  EXPECT_NEAR(getter::element(m33, 1u, 0u), 6.f, this->tolerance());
  EXPECT_NEAR(getter::element(m33, 1u, 1u), 8.f, this->tolerance());
  EXPECT_NEAR(getter::element(m33, 1u, 2u), 10.f, this->tolerance());
  EXPECT_NEAR(getter::element(m33, 2u, 0u), 9.f, this->tolerance());
  EXPECT_NEAR(getter::element(m33, 2u, 1u), 12.f, this->tolerance());
  EXPECT_NEAR(getter::element(m33, 2u, 2u), 15.f, this->tolerance());
}

TEST_F(detray_algebra, cholesky_decomposition) {
  using matrix3_type = matrix<3, 3>;

  // Define A
  auto A = detray::matrix::zero<matrix3_type>();
  getter::element(A, 0u, 0u) = 4.f;
  getter::element(A, 0u, 1u) = 12.f;
  getter::element(A, 0u, 2u) = -16.f;
  getter::element(A, 1u, 0u) = 12.f;
  getter::element(A, 1u, 1u) = 37.f;
  getter::element(A, 1u, 2u) = -43.f;
  getter::element(A, 2u, 0u) = -16.f;
  getter::element(A, 2u, 1u) = -43.f;
  getter::element(A, 2u, 2u) = 98.f;

  // Get L that satisfies A = L * L^T and check if it is the expected value
  const matrix3_type L = detray::matrix::cholesky_decomposition(A);

  EXPECT_FLOAT_EQ(static_cast<float>(getter::element(L, 0u, 0u)), 2.f);
  EXPECT_FLOAT_EQ(static_cast<float>(getter::element(L, 0u, 1u)), 0.f);
  EXPECT_FLOAT_EQ(static_cast<float>(getter::element(L, 0u, 2u)), 0.f);
  EXPECT_FLOAT_EQ(static_cast<float>(getter::element(L, 1u, 0u)), 6.f);
  EXPECT_FLOAT_EQ(static_cast<float>(getter::element(L, 1u, 1u)), 1.f);
  EXPECT_FLOAT_EQ(static_cast<float>(getter::element(L, 1u, 2u)), 0.f);
  EXPECT_FLOAT_EQ(static_cast<float>(getter::element(L, 2u, 0u)), -8.f);
  EXPECT_FLOAT_EQ(static_cast<float>(getter::element(L, 2u, 1u)), 5.f);
  EXPECT_FLOAT_EQ(static_cast<float>(getter::element(L, 2u, 2u)), 3.f);

  // Compare A and L * L^T
  const matrix3_type B = L * detray::matrix::transpose(L);

  for (unsigned int i = 0u; i < 3u; i++) {
    for (unsigned int j = 0u; j < 3u; j++) {
      EXPECT_FLOAT_EQ(static_cast<float>(getter::element(A, i, j)),
                      static_cast<float>(getter::element(B, i, j)));
    }
  }
}

}  // namespace detray::test
