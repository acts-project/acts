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
#include "algebra/fastor.hpp"
#include "detray/algebra/concepts.hpp"

// System include(s)
#include <algorithm>
#include <vector>

namespace detray::algebra {

/// Fill an @c Fastor based vector with random values
template <detray::concepts::vector vector_t>
inline void fill_random_vec(std::vector<vector_t>& collection) {
  auto rand_obj = [](vector_t& v) { v.random(); };

  collection.resize(collection.capacity());
  std::ranges::for_each(collection, rand_obj);
}

/// Fill a @c Fastor based transform3 with random values
template <detray::concepts::transform3D transform3_t>
inline void fill_random_trf(std::vector<transform3_t>& collection) {
  using vector_t = typename transform3_t::vector3;

  auto rand_obj = [](transform3_t& trf) {
    vector_t x_axis;
    vector_t z_axis;
    vector_t t;

    x_axis.random();
    x_axis = detray::vector::normalize(x_axis);
    z_axis.random();
    t.random();
    t = detray::vector::normalize(t);

    // Gram-Schmidt projection
    typename transform3_t::scalar_type coeff =
        detray::vector::dot(x_axis, z_axis) / detray::vector::norm(x_axis);
    z_axis = x_axis - coeff * z_axis;

    trf = transform3_t{t, x_axis, detray::vector::normalize(z_axis)};
  };

  collection.resize(collection.capacity());
  std::ranges::for_each(collection, rand_obj);
}

/// Fill a @c Fastor based matrix with random values
template <detray::concepts::matrix matrix_t>
inline void fill_random_matrix(std::vector<matrix_t>& collection) {
  auto rand_obj = [](matrix_t& m) { m.random(); };

  collection.resize(collection.capacity());
  std::ranges::for_each(collection, rand_obj);
}

}  // namespace detray::algebra
