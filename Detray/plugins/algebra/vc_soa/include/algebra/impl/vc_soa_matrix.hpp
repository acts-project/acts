// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "algebra/impl/vc_soa_vector.hpp"
#include "detray/algebra/common/matrix.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

namespace detray::algebra::vc_soa::math {

using algebra::storage::identity;
using algebra::storage::set_identity;
using algebra::storage::set_zero;
using algebra::storage::transpose;
using algebra::storage::zero;

template <std::size_t ROW, std::size_t COL, concepts::simd_scalar scalar_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr scalar_t determinant(
    const algebra::storage::matrix<array_t, scalar_t, ROW, COL>
        & /*unused*/) noexcept {
  // @TODO: Implement
  return static_cast<scalar_t>(0);
}

template <std::size_t ROW, std::size_t COL, concepts::simd_scalar scalar_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr algebra::storage::matrix<array_t, scalar_t, ROW,
                                                      COL>
inverse(
    const algebra::storage::matrix<array_t, scalar_t, ROW, COL> &m) noexcept {
  // @TODO: Implement
  return m;
}

}  // namespace detray::algebra::vc_soa::math
