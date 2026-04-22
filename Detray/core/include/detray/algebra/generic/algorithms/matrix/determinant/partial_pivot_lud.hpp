// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/algebra/concepts.hpp"
#include "detray/algebra/generic/algorithms/matrix/decomposition/partial_pivot_lud.hpp"
#include "detray/algebra/type_traits.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

namespace detray::algebra::generic::matrix::determinant {

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

  DETRAY_HOST_DEVICE constexpr scalar_t operator()(const matrix_t& m) const {
    constexpr element_getter_t elem{};
    constexpr index_t N{detray::traits::max_rank<matrix_t>};

    const typename decomposition_t::template lud<
        detray::traits::max_rank<matrix_t>>
        decomp_res = decomposition_t()(m);

    // Get the LU decomposition matrix equal to (L - I) + U
    const auto& lu = decomp_res.lu;
    const auto n_pivot = static_cast<index_t>(decomp_res.n_pivot);

    scalar_t det = elem(lu, 0, 0);

    for (index_t i = 1; i < N; i++) {
      det *= elem(lu, i, i);
    }

    return (n_pivot - N) % 2 == 0 ? det : -det;
  }
};

}  // namespace detray::algebra::generic::matrix::determinant
