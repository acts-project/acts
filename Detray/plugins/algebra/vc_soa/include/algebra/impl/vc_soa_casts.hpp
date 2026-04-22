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

// Vc include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#include <Vc/Vc>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC

namespace detray::algebra {

// Forward declare the generic cast impl from matrices
template <concepts::value value_t, concepts::matrix matrix_t>
DETRAY_HOST_DEVICE constexpr auto cast_to(const matrix_t& m);

/// Cast a Vc SoA salar @param s to the precision given by @tparam other_value_t
template <concepts::value other_value_t, concepts::value value_t>
DETRAY_HOST constexpr auto cast_to(const Vc::Vector<value_t>& s) {
  using other_scalar_t = Vc::Vector<other_value_t>;

  return Vc::simd_cast<other_scalar_t>(s);
}

/// Cast a Vc SoA transform @param trf to the precision given by @tparam value_t
template <concepts::value value_t, concepts::transform3D transform_t>
  requires Vc::is_simd_vector<typename transform_t::scalar_type>::value
DETRAY_HOST_DEVICE constexpr auto cast_to(const transform_t& trf) {
  using scalar_t = Vc::Vector<value_t>;
  using new_trf3_t = typename transform_t::template other_type<scalar_t>;

  return new_trf3_t{cast_to<value_t>(trf.matrix()),
                    cast_to<value_t>(trf.matrix_inverse())};
}

}  // namespace detray::algebra
