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
#include "detray/algebra/type_traits.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s)
#include <type_traits>

namespace detray::algebra::generic::matrix::determinant {

/// Evaluate determinants with direct formulas for supported matrix ranks.
template <concepts::square_matrix matrix_t>
struct hard_coded {
  /// Matrix element type.
  using scalar_t = detray::traits::value_t<matrix_t>;
  /// Matrix index type.
  using index_t = detray::traits::index_t<matrix_t>;

  /// Function (object) used for accessing a matrix element
  using element_getter_t = detray::traits::element_getter_t<matrix_t>;

  /// Evaluate the determinant of a 2 by 2 matrix.
  /// @param m Matrix whose determinant is calculated.
  /// @returns The determinant.
  template <typename M = matrix_t>
    requires(detray::traits::max_rank<M> == 2)
  DETRAY_HOST_DEVICE constexpr scalar_t operator()(const matrix_t &m) const {
    constexpr element_getter_t elem{};

    return elem(m, 0, 0) * elem(m, 1, 1) - elem(m, 0, 1) * elem(m, 1, 0);
  }

  /// Evaluate the determinant of a 4 by 4 matrix.
  /// @param m Matrix whose determinant is calculated.
  /// @returns The determinant.
  template <typename M = matrix_t>
    requires(detray::traits::max_rank<M> == 4)
  DETRAY_HOST_DEVICE constexpr scalar_t operator()(const matrix_t &m) const {
    constexpr element_getter_t elem{};

    return elem(m, 0, 3) * elem(m, 1, 2) * elem(m, 2, 1) * elem(m, 3, 0) -
           elem(m, 0, 2) * elem(m, 1, 3) * elem(m, 2, 1) * elem(m, 3, 0) -
           elem(m, 0, 3) * elem(m, 1, 1) * elem(m, 2, 2) * elem(m, 3, 0) +
           elem(m, 0, 1) * elem(m, 1, 3) * elem(m, 2, 2) * elem(m, 3, 0) +
           elem(m, 0, 2) * elem(m, 1, 1) * elem(m, 2, 3) * elem(m, 3, 0) -
           elem(m, 0, 1) * elem(m, 1, 2) * elem(m, 2, 3) * elem(m, 3, 0) -
           elem(m, 0, 3) * elem(m, 1, 2) * elem(m, 2, 0) * elem(m, 3, 1) +
           elem(m, 0, 2) * elem(m, 1, 3) * elem(m, 2, 0) * elem(m, 3, 1) +
           elem(m, 0, 3) * elem(m, 1, 0) * elem(m, 2, 2) * elem(m, 3, 1) -
           elem(m, 0, 0) * elem(m, 1, 3) * elem(m, 2, 2) * elem(m, 3, 1) -
           elem(m, 0, 2) * elem(m, 1, 0) * elem(m, 2, 3) * elem(m, 3, 1) +
           elem(m, 0, 0) * elem(m, 1, 2) * elem(m, 2, 3) * elem(m, 3, 1) +
           elem(m, 0, 3) * elem(m, 1, 1) * elem(m, 2, 0) * elem(m, 3, 2) -
           elem(m, 0, 1) * elem(m, 1, 3) * elem(m, 2, 0) * elem(m, 3, 2) -
           elem(m, 0, 3) * elem(m, 1, 0) * elem(m, 2, 1) * elem(m, 3, 2) +
           elem(m, 0, 0) * elem(m, 1, 3) * elem(m, 2, 1) * elem(m, 3, 2) +
           elem(m, 0, 1) * elem(m, 1, 0) * elem(m, 2, 3) * elem(m, 3, 2) -
           elem(m, 0, 0) * elem(m, 1, 1) * elem(m, 2, 3) * elem(m, 3, 2) -
           elem(m, 0, 2) * elem(m, 1, 1) * elem(m, 2, 0) * elem(m, 3, 3) +
           elem(m, 0, 1) * elem(m, 1, 2) * elem(m, 2, 0) * elem(m, 3, 3) +
           elem(m, 0, 2) * elem(m, 1, 0) * elem(m, 2, 1) * elem(m, 3, 3) -
           elem(m, 0, 0) * elem(m, 1, 2) * elem(m, 2, 1) * elem(m, 3, 3) -
           elem(m, 0, 1) * elem(m, 1, 0) * elem(m, 2, 2) * elem(m, 3, 3) +
           elem(m, 0, 0) * elem(m, 1, 1) * elem(m, 2, 2) * elem(m, 3, 3);
  }
};

}  // namespace detray::algebra::generic::matrix::determinant
