// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/algebra/concepts.hpp"
#include "detray/algebra/generic/algorithms/matrix/decomposition/partial_pivot_lud.hpp"
#include "detray/algebra/type_traits.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

namespace detray::algebra::generic::matrix::inverse {

/// "Partial Pivot LU Decomposition", assuming a N X N matrix
template <concepts::square_matrix matrix_t>
struct partial_pivot_lud {
  using scalar_t = detray::traits::value_t<matrix_t>;
  using index_t = detray::traits::index_t<matrix_t>;

  /// Function (object) used for accessing a matrix element
  using element_getter_t = detray::traits::element_getter_t<matrix_t>;

  using decomposition_t =
      typename algebra::generic::matrix::decomposition::partial_pivot_lud<
          matrix_t>;

  DETRAY_HOST_DEVICE constexpr matrix_t operator()(const matrix_t& m) const {
    constexpr element_getter_t elem{};
    constexpr index_t N{detray::traits::max_rank<matrix_t>};

    const typename decomposition_t::template lud<N> decomp_res =
        decomposition_t()(m);

    // Get the LU decomposition matrix equal to (L - I) + U
    const auto& lu = decomp_res.lu;

    // Permutation vector
    const auto& P = decomp_res.P;

    // Inverse matrix
    matrix_t inv{};

    // Calculate inv(A) = inv(U) * inv(L) * P;
    for (index_t j = 0; j < N; j++) {
      for (index_t i = 0; i < N; i++) {
        elem(inv, i, j) = static_cast<index_t>(P[i]) == j
                              ? static_cast<scalar_t>(1.0)
                              : static_cast<scalar_t>(0.0);

        for (index_t k = 0; k < i; k++) {
          elem(inv, i, j) -= elem(lu, i, k) * elem(inv, k, j);
        }
      }

      for (index_t i = N - 1; int(i) >= 0; i--) {
        for (index_t k = i + 1; k < N; k++) {
          elem(inv, i, j) -= elem(lu, i, k) * elem(inv, k, j);
        }
        elem(inv, i, j) /= elem(lu, i, i);
      }
    }

    return inv;
  }
};

}  // namespace detray::algebra::generic::matrix::inverse
