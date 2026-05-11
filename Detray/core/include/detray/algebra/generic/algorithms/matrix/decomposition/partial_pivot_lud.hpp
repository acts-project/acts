// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/algebra/common/math.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s).
#include <array>

namespace detray::algebra::generic::matrix::decomposition {

/// "Partial Pivot LU Decomposition", assuming a N X N matrix
template <concepts::matrix matrix_t>
struct partial_pivot_lud {
  using scalar_t = detray::traits::value_t<matrix_t>;
  using index_t = detray::traits::index_t<matrix_t>;
  using vector_t = detray::traits::vector_t<matrix_t>;

  template <index_t N>
  struct lud {
    // LU decomposition matrix, equal to (L - I) + U, where the diagonal
    // components of L is always 1
    matrix_t lu;

    // Permutation vector
    vector_t P;

    // Number of pivots
    int n_pivot = 0;
  };

  DETRAY_HOST_DEVICE constexpr lud<detray::traits::max_rank<matrix_t>>
  operator()(const matrix_t& m) const {
    // Function (object) used for accessing a matrix element
    using element_getter_t = detray::traits::element_getter_t<matrix_t>;

    constexpr element_getter_t elem{};
    constexpr index_t N{detray::traits::max_rank<matrix_t>};

    // LU decomposition matrix
    matrix_t lu = m;

    // Permutation
    vector_t P;

    // Max index and value
    index_t max_idx;
    scalar_t max_val;
    scalar_t abs_val;

    // Number of pivoting
    int n_pivot = N;

    // Rows for swapping
    vector_t row_0;
    vector_t row_1;

    // Unit permutation matrix, P[N] initialized with N
    for (index_t i = 0; i < N; i++) {
      P[i] = static_cast<scalar_t>(i);
    }

    for (index_t i = 0; i < N; i++) {
      max_val = 0;
      max_idx = i;

      for (index_t k = i; k < N; k++) {
        abs_val = algebra::math::fabs(elem(lu, k, i));

        if (abs_val > max_val) {
          max_val = abs_val;
          max_idx = k;
        }
      }

      if (max_idx != i) {
        // Pivoting P
        auto j = P[i];

        P[i] = P[max_idx];
        P[max_idx] = j;

        // Pivoting rows of A
        for (index_t q = 0; q < N; q++) {
          row_0[q] = elem(lu, i, q);
          row_1[q] = elem(lu, max_idx, q);
        }
        for (index_t q = 0; q < N; q++) {
          elem(lu, i, q) = row_1[q];
          elem(lu, max_idx, q) = row_0[q];
        }

        // counting pivots starting from N (for determinant)
        n_pivot++;
      }

      for (index_t j = i + 1; j < N; j++) {
        // m[j][i] /= m[i][i];
        elem(lu, j, i) /= elem(lu, i, i);

        for (index_t k = i + 1; k < N; k++) {
          // m[j][k] -= m[j][i] * m[i][k];
          elem(lu, j, k) -= elem(lu, j, i) * elem(lu, i, k);
        }
      }
    }

    return {lu, P, n_pivot};
  }
};

}  // namespace detray::algebra::generic::matrix::decomposition
