// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/algebra/common/matrix.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/algebra/generic/generic.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

namespace detray::algebra::vc_aos::math {

using algebra::storage::identity;
using algebra::storage::set_identity;
using algebra::storage::set_zero;
using algebra::storage::transpose;
using algebra::storage::zero;

/// @returns the determinant
template <std::size_t N, concepts::value value_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr value_t determinant(
    const algebra::storage::matrix<array_t, value_t, N, N> &m) noexcept {
  // TODO: Implement vectorization friendly version
  return algebra::generic::math::determinant(m);
}

/// @returns the inverse
template <std::size_t ROW, std::size_t COL, concepts::value value_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr algebra::storage::matrix<array_t, value_t, COL,
                                                      ROW>
inverse(
    const algebra::storage::matrix<array_t, value_t, ROW, COL> &m) noexcept {
  // TODO: Implement vectorization friendly version
  return algebra::generic::math::inverse(m);
}

/// @returns the transpose
template <std::size_t ROW, std::size_t COL, concepts::value value_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr algebra::storage::matrix<array_t, value_t, COL,
                                                      ROW>
transpose(
    const algebra::storage::matrix<array_t, value_t, ROW, COL> &m) noexcept {
  // TODO: Implement vectorization friendly version
  return algebra::generic::math::transpose(m);
}

}  // namespace detray::algebra::vc_aos::math
