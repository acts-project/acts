// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/algebra/common/matrix.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s).
#include <cstddef>

/// Matrix functions for the SoA plugins (Vc SoA and array SoA)
namespace detray::algebra::soa {

using algebra::storage::identity;
using algebra::storage::set_identity;
using algebra::storage::set_zero;
using algebra::storage::transpose;
using algebra::storage::zero;

namespace detail {

/// Dependent false value to trigger a static assertion on instantiation
template <typename>
inline constexpr bool always_false_v = false;

}  // namespace detail

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
    const algebra::storage::matrix<array_t, scalar_t, ROW,
                                   COL>& /*m*/) noexcept {
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
inverse(const algebra::storage::matrix<array_t, scalar_t, ROW,
                                       COL>& /*m*/) noexcept {
  static_assert(detail::always_false_v<scalar_t>,
                "inverse is not implemented for the SoA plugins: the generic "
                "algorithms return traits::value_t, not the SIMD type");
  return {};
}

}  // namespace detray::algebra::soa
