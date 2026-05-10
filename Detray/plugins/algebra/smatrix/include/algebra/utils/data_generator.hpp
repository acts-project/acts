// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "algebra/smatrix.hpp"
#include "detray/algebra/concepts.hpp"

// System include(s)
#include <algorithm>
#include <random>
#include <vector>

namespace detray::algebra {

/// Fill an @c std::array based vector with random values
template <detray::concepts::vector vector_t>
inline void fill_random_vec(std::vector<vector_t> &collection) {
  // Generate a vector of the right type with random values
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<typename vector_t::value_type> dist(0.f, 1.f);

  auto rand_obj = [&]() {
    vector_t tmp{};

    for (unsigned int i = 0u; i < detray::traits::size<vector_t>; ++i) {
      tmp[i] = dist(mt);
    }

    return tmp;
  };

  collection.resize(collection.capacity());
  std::ranges::generate(collection, rand_obj);
}

/// Fill a @c std::array based transform3 with random values
template <detray::concepts::transform3D transform3_t>
inline void fill_random_trf(std::vector<transform3_t> &collection) {
  using vector_t = typename transform3_t::vector3;

  // Generate a random, but valid affine transformation
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<typename transform3_t::scalar_type> dist(0.f,
                                                                          1.f);

  auto rand_obj = [&]() {
    vector_t x_axis;
    vector_t z_axis;
    vector_t t;

    x_axis = detray::vector::normalize(vector_t{dist(mt), dist(mt), dist(mt)});
    z_axis = {dist(mt), dist(mt), dist(mt)};
    t = detray::vector::normalize(vector_t{dist(mt), dist(mt), dist(mt)});

    // Gram-Schmidt projection
    typename transform3_t::scalar_type coeff =
        detray::vector::dot(x_axis, z_axis) / detray::vector::norm(x_axis);
    z_axis = x_axis - coeff * z_axis;

    return transform3_t{t, x_axis, detray::vector::normalize(z_axis)};
  };

  collection.resize(collection.capacity());
  std::ranges::generate(collection, rand_obj);
}

/// Fill a @c std::array based matrix with random values
template <detray::concepts::matrix matrix_t>
inline void fill_random_matrix(std::vector<matrix_t> &collection) {
  using scalar_t = typename matrix_t::value_type;

  // Generate a random, but valid affine transformation
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<scalar_t> dist(0.f, 1.f);
  auto rand_obj = [&]() {
    matrix_t m;

    for (unsigned int j = 0u; j < detray::traits::columns<matrix_t>; ++j) {
      for (unsigned int i = 0u; i < detray::traits::rows<matrix_t>; ++i) {
        m(i, j) = dist(mt);
      }
    }

    return m;
  };

  collection.resize(collection.capacity());
  std::ranges::generate(collection, rand_obj);
}

}  // namespace detray::algebra
