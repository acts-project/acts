// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "algebra/vc_soa.hpp"
#include "detray/algebra/concepts.hpp"

// System include(s)
#include <algorithm>
#include <vector>

namespace detray::algebra {

/// Fill a @c Vc::Vector based vector with random values
template <detray::concepts::vector vector_soa_t>
inline void fill_random_vec(std::vector<vector_soa_t> &collection) {
  // Generate a vector of the right type with random values
  auto rand_obj = []() {
    using simd_vector_t = typename vector_soa_t::scalar_type;
    vector_soa_t tmp{};

    for (std::size_t i = 0u; i < detray::traits::size<vector_soa_t>; ++i) {
      tmp[i] = simd_vector_t::Random();
    }

    return tmp;
  };

  collection.resize(collection.capacity());
  std::ranges::generate(collection, rand_obj);
}

/// Fill a @c Vc::Vector based transform3 with random values
template <detray::concepts::transform3D transform3_t>
inline void fill_random_trf(std::vector<transform3_t> &collection) {
  // Generate a random, but valid affine transformation
  auto rand_obj = []() {
    using vector_t = typename transform3_t::vector3;
    using simd_vector_t = typename transform3_t::scalar_type;

    vector_t x_axis;
    vector_t z_axis;
    vector_t t;

    x_axis[0] = simd_vector_t::Random();
    x_axis[1] = simd_vector_t::Random();
    x_axis[2] = simd_vector_t::Random();
    x_axis = detray::vector::normalize(x_axis);

    z_axis[0] = simd_vector_t::Random();
    z_axis[1] = simd_vector_t::Random();
    z_axis[2] = simd_vector_t::Random();

    t[0] = simd_vector_t::Random();
    t[1] = simd_vector_t::Random();
    t[2] = simd_vector_t::Random();
    t = detray::vector::normalize(t);

    // Gram-Schmidt projection
    simd_vector_t coeff =
        detray::vector::dot(x_axis, z_axis) / detray::vector::norm(x_axis);
    z_axis = x_axis - coeff * z_axis;

    return transform3_t{t, x_axis, detray::vector::normalize(z_axis)};
  };

  collection.resize(collection.capacity());
  std::ranges::generate(collection, rand_obj);
}

/// Fill a @c Vc::Vector based matrix with random values
template <detray::concepts::matrix matrix_t>
inline void fill_random_matrix(std::vector<matrix_t> &collection) {
  // Generate a random, but valid affine transformation
  auto rand_obj = []() {
    using simd_vector_t = typename matrix_t::scalar_type;

    matrix_t m;

    for (std::size_t j = 0u; j < matrix_t::columns(); ++j) {
      typename matrix_t::vector_type v;

      for (std::size_t i = 0u; i < matrix_t::rows(); ++i) {
        v[i] = simd_vector_t::Random();
      }

      m[j] = v;
    }

    return m;
  };

  collection.resize(collection.capacity());
  std::ranges::generate(collection, rand_obj);
}

}  // namespace detray::algebra
