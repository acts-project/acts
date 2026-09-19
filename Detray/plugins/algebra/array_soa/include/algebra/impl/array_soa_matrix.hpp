// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "algebra/impl/array_soa_vector.hpp"
#include "detray/algebra/common/matrix.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

namespace detray::algebra::array_soa::math {

using algebra::storage::identity;
using algebra::storage::set_identity;
using algebra::storage::set_zero;
using algebra::storage::transpose;
using algebra::storage::zero;

namespace detail {

template <typename>
inline constexpr bool always_false_v = false;

}  // namespace detail

template <std::size_t ROW, std::size_t COL, concepts::array_soa_simd scalar_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr scalar_t determinant(
    const algebra::storage::matrix<array_t, scalar_t, ROW,
                                   COL>& /*m*/) noexcept {
  static_assert(
      detail::always_false_v<scalar_t>,
      "determinant is not implemented for the SoA plugin: the generic "
      "algorithms return traits::value_t, not the lane bundle");
  return {};
}

template <std::size_t ROW, std::size_t COL, concepts::array_soa_simd scalar_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr algebra::storage::matrix<array_t, scalar_t, ROW,
                                                      COL>
inverse(const algebra::storage::matrix<array_t, scalar_t, ROW,
                                       COL>& /*m*/) noexcept {
  static_assert(detail::always_false_v<scalar_t>,
                "inverse is not implemented for the SoA plugin: the generic "
                "algorithms return traits::value_t, not the lane bundle");
  return {};
}

}  // namespace detray::algebra::array_soa::math
