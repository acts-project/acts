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
#include "detray/algebra/common/math.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/algebra/generic/generic.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

namespace detray::algebra::array {

/// @returns zero matrix of type @tparam matrix_t
template <concepts::matrix matrix_t>
  requires(std::is_scalar_v<typename matrix_t::value_type::value_type>)
DETRAY_HOST_DEVICE constexpr matrix_t zero() {
  matrix_t ret;

  for (std::size_t j = 0; j < detray::traits::columns<matrix_t>; ++j) {
    for (std::size_t i = 0; i < detray::traits::rows<matrix_t>; ++i) {
      ret[j][i] = 0;
    }
  }

  return ret;
}

/// @returns identity matrix of type @tparam matrix_t
template <concepts::matrix matrix_t>
  requires(std::is_scalar_v<typename matrix_t::value_type::value_type>)
DETRAY_HOST_DEVICE constexpr auto identity() {
  auto ret{zero<matrix_t>()};

  for (std::size_t i = 0; i < detray::traits::max_rank<matrix_t>; ++i) {
    ret[i][i] = 1;
  }

  return ret;
}

/// Set @param m as zero matrix
template <std::size_t ROWS, std::size_t COLS, concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr void set_zero(
    array_t<array_t<scalar_t, ROWS>, COLS> &m) {
  m = zero<array_t<array_t<scalar_t, ROWS>, COLS>>();
}

/// Set @param m as identity matrix
template <std::size_t ROWS, std::size_t COLS, concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr void set_identity(
    array_t<array_t<scalar_t, ROWS>, COLS> &m) {
  m = identity<array_t<array_t<scalar_t, ROWS>, COLS>>();
}

/// @returns the transpose matrix of @param m
template <std::size_t ROWS, std::size_t COLS, concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr auto transpose(
    const array_t<array_t<scalar_t, ROWS>, COLS> &m) {
  return algebra::generic::math::transpose(m);
}

/// @returns the determinant of @param m
template <std::size_t ROWS, std::size_t COLS, concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr scalar_t determinant(
    const array_t<array_t<scalar_t, ROWS>, COLS> &m) {
  return algebra::generic::math::determinant(m);
}

/// @returns the determinant of @param m
template <std::size_t ROWS, std::size_t COLS, concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr auto inverse(
    const array_t<array_t<scalar_t, ROWS>, COLS> &m) {
  return algebra::generic::math::inverse(m);
}

}  // namespace detray::algebra::array
