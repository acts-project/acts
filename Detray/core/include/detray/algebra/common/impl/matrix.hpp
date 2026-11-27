// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/algebra/common/storage/matrix.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s).
#include <cstddef>

// Matrix functions for the AoS/SoA plugins
//
// This header only blocks the generic determinant and inverse, which are not
// correct for SIMD scalars.
namespace detray::algebra::common {

namespace detail {

/// Dependent false value to trigger a static assertion on instantiation
template <typename>
inline constexpr bool always_false_v = false;

}  // namespace detail

/// Get a zero-initialized matrix
template <concepts::matrix matrix_t, std::size_t COLS = matrix_t::columns()>
DETRAY_HOST_DEVICE constexpr matrix_t zero() noexcept {
  matrix_t m;

  DETRAY_UNROLL_N(COLS)
  for (std::size_t j = 0u; j < COLS; ++j) {
    // Fill zero initialized vector
    m[j] = typename matrix_t::vector_type{};
  }

  return m;
}

/// Set a matrix to zero
template <concepts::matrix matrix_t>
DETRAY_HOST_DEVICE constexpr void set_zero(matrix_t &m) noexcept {
  m = zero<matrix_t>();
}

/// Build an identity matrix
template <concepts::matrix matrix_t,
          std::size_t R = detray::traits::max_rank<matrix_t>>
DETRAY_HOST_DEVICE constexpr matrix_t identity() noexcept {
  // Zero initialized
  matrix_t m{zero<matrix_t>()};

  DETRAY_UNROLL_N(R)
  for (std::size_t i = 0u; i < R; ++i) {
    m[i][i] = typename matrix_t::scalar_type(1);
  }

  return m;
}

/// Set a matrix to zero
template <concepts::matrix matrix_t>
DETRAY_HOST_DEVICE constexpr void set_identity(matrix_t &m) noexcept {
  m = identity<matrix_t>();
}

/// Transpose the matrix @param m
template <std::size_t ROW, std::size_t COL, concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t, std::size_t... I>
DETRAY_HOST_DEVICE constexpr auto transpose(
    const storage::matrix<array_t, scalar_t, ROW, COL> &m,
    std::index_sequence<I...> /*unused*/) noexcept {
  using matrix_T_t = storage::matrix<array_t, scalar_t, COL, ROW>;
  using column_t = typename matrix_T_t::vector_type;

  matrix_T_t res_m;

  DETRAY_UNROLL_N(ROW)
  for (std::size_t j = 0u; j < ROW; ++j) {
    res_m[j] = column_t{m[I][j]...};
  }

  return res_m;
}

/// Build an identity matrix
template <concepts::matrix matrix_t>
DETRAY_HOST_DEVICE constexpr auto transpose(const matrix_t &m) noexcept {
  return transpose(m, std::make_index_sequence<matrix_t::columns()>());
}

/// Determinant of a matrix with SIMD elements: not implemented
///
/// @note The generic algorithms return @c traits::value_t instead of the SIMD
/// scalar type. The function fails at compile time, so that a wrong result is
/// never returned.
///
/// @tparam ROW Number of matrix rows
/// @tparam COL Number of matrix columns
/// @tparam scalar_t SIMD scalar type
/// @tparam array_t Underlying storage array template
///
/// @returns no value, the instantiation fails
template <std::size_t ROW, std::size_t COL, concepts::simd_scalar scalar_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr scalar_t determinant(
    const algebra::storage::matrix<array_t, scalar_t, ROW, COL>
        & /*m*/) noexcept {
  static_assert(detail::always_false_v<scalar_t>,
                "determinant is not implemented for the SoA plugins: the "
                "generic algorithms return traits::value_t, not the SIMD type");
  return {};
}

/// Inverse of a matrix with SIMD elements: not implemented
///
/// @note The generic algorithms return @c traits::value_t instead of the SIMD
/// scalar type. The function fails at compile time, so that a wrong result is
/// never returned.
///
/// @tparam ROW Number of matrix rows
/// @tparam COL Number of matrix columns
/// @tparam scalar_t SIMD scalar type
/// @tparam array_t Underlying storage array template
///
/// @returns no value, the instantiation fails
template <std::size_t ROW, std::size_t COL, concepts::simd_scalar scalar_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr algebra::storage::matrix<array_t, scalar_t, ROW,
                                                      COL>
inverse(const algebra::storage::matrix<array_t, scalar_t, ROW, COL>
            & /*m*/) noexcept {
  static_assert(detail::always_false_v<scalar_t>,
                "inverse is not implemented for the SoA plugins: the generic "
                "algorithms return traits::value_t, not the SIMD type");
  return {};
}

}  // namespace detray::algebra::common
