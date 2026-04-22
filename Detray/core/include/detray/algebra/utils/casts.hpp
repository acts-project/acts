// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/algebra/concepts.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s)
#include <concepts>

namespace detray::algebra {

/// Cast a salar (might be simd) @param s to the precision given by
/// @tparam value_t
template <concepts::value value_t, concepts::scalar scalar_t>
  requires std::convertible_to<scalar_t, value_t>
DETRAY_HOST_DEVICE constexpr auto cast_to(const scalar_t& s) {
  return static_cast<value_t>(s);
}

/// Cast a generic vector or point @param v to the precision given by
/// @tparam value_t
template <concepts::value value_t, typename vector_t>
  requires(concepts::vector<vector_t> || concepts::point<vector_t>)
DETRAY_HOST_DEVICE constexpr auto cast_to(const vector_t& v) {
  using index_t = detray::traits::index_t<vector_t>;

  constexpr index_t size{detray::traits::size<vector_t>};

  using new_vector_t = detray::traits::get_vector_t<vector_t, size, value_t>;
  new_vector_t ret;

  static_assert(std::same_as<value_t, detray::traits::value_t<new_vector_t>>);

  for (index_t i = 0; i < size; ++i) {
    ret[i] = ::detray::algebra::cast_to<value_t>(v[i]);
  }

  return ret;
}

/// Cast a column matrix @param v to the precision given by @tparam value_t
template <concepts::value value_t, concepts::column_matrix vector_t>
DETRAY_HOST_DEVICE constexpr auto cast_to(const vector_t& v) {
  using index_t = detray::traits::index_t<vector_t>;
  using element_getter_t = detray::traits::element_getter_t<vector_t>;

  constexpr index_t rows{detray::traits::rows<vector_t>};

  using new_vector_t = detray::traits::get_matrix_t<vector_t, rows, 1, value_t>;
  new_vector_t ret;

  static_assert(std::same_as<value_t, detray::traits::value_t<new_vector_t>>);

  for (index_t i = 0; i < rows; ++i) {
    element_getter_t{}(ret, i, 0) =
        ::detray::algebra::cast_to<value_t>(element_getter_t{}(v, i, 0));
  }

  return ret;
}

/// Cast a generic matrix @param v to the precision given by @tparam value_t
template <concepts::value value_t, concepts::matrix matrix_t>
DETRAY_HOST_DEVICE constexpr auto cast_to(const matrix_t& m) {
  using index_t = detray::traits::index_t<matrix_t>;
  using element_getter_t = detray::traits::element_getter_t<matrix_t>;

  constexpr index_t rows{detray::traits::rows<matrix_t>};
  constexpr index_t columns{detray::traits::columns<matrix_t>};

  using new_matrix_t =
      detray::traits::get_matrix_t<matrix_t, rows, columns, value_t>;
  new_matrix_t ret;

  static_assert(std::same_as<value_t, detray::traits::value_t<new_matrix_t>>);

  for (index_t j = 0; j < columns; ++j) {
    for (index_t i = 0; i < rows; ++i) {
      element_getter_t{}(ret, i, j) =
          ::detray::algebra::cast_to<value_t>(element_getter_t{}(m, i, j));
    }
  }

  return ret;
}

/// Cast a 3D transform @param trf to the precision given by @tparam scalar_t
template <concepts::scalar scalar_t, concepts::transform3D transform_t>
  requires(!concepts::simd_scalar<scalar_t> &&
           !concepts::simd_scalar<typename transform_t::scalar_type>)
DETRAY_HOST_DEVICE constexpr auto cast_to(const transform_t& trf) {
  using new_trf3_t = typename transform_t::template other_type<scalar_t>;

  return new_trf3_t{::detray::algebra::cast_to<scalar_t>(trf.matrix()),
                    ::detray::algebra::cast_to<scalar_t>(trf.matrix_inverse())};
}

}  // namespace detray::algebra
