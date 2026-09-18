// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "algebra/impl/array_soa_simd.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s)
#include <concepts>
#include <cstddef>

namespace detray::algebra {

// Forward declare the generic cast impl from matrices
template <concepts::value value_t, concepts::matrix matrix_t>
DETRAY_HOST_DEVICE constexpr auto cast_to(const matrix_t& m);

/// Cast a lane bundle @param s to the precision given by @tparam other_value_t
template <concepts::value other_value_t, concepts::value value_t, std::size_t W>
DETRAY_HOST_DEVICE constexpr auto cast_to(
    const array_soa::simd<value_t, W>& s) {
  if constexpr (std::same_as<other_value_t, value_t>) {
    return s;
  } else {
    return array_soa::simd<other_value_t, W>{s};
  }
}

/// Cast a lane bundle based transform @param trf to the precision given by @tparam value_t
template <concepts::value value_t, concepts::transform3D transform_t>
  requires concepts::array_soa_simd<typename transform_t::scalar_type>
DETRAY_HOST_DEVICE constexpr auto cast_to(const transform_t& trf) {
  using scalar_t = array_soa::simd<value_t, transform_t::scalar_type::size()>;
  using new_trf3_t = typename transform_t::template other_type<scalar_t>;

  return new_trf3_t{cast_to<value_t>(trf.matrix()),
                    cast_to<value_t>(trf.matrix_inverse())};
}

}  // namespace detray::algebra
