// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "algebra/array_soa.hpp"
#include "detray/algebra/concepts.hpp"

// System include(s)
#include <algorithm>
#include <random>
#include <vector>

namespace detray::algebra {

/// Fill an array SoA with random values
template <detray::concepts::scalar scalar_t>
inline void fill_random_scalar(std::vector<scalar_t> &collection) {
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(0.f, 1.f);

  auto rand_simd = [&]() {
    scalar_t s;
    for (std::size_t i = 0u; i < scalar_t::size(); ++i) {
      s[i] = dist(mt);
    }
    return static_cast<scalar_t>(s);
  };

  collection.resize(collection.capacity());
  std::ranges::generate(collection, rand_obj);
}

/// Fill a lane bundle based vector with random values
template <detray::concepts::vector vector_soa_t>
inline void fill_random_vec(std::vector<vector_soa_t> &collection) {
  using simd_vector_t = typename vector_soa_t::scalar_type;
  using value_t = typename simd_vector_t::value_type;

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<value_t> dist(0.f, 1.f);

  auto rand_simd = [&]() {
    simd_vector_t s;
    for (std::size_t i = 0u; i < simd_vector_t::size(); ++i) {
      s[i] = dist(mt);
    }
    return s;
  };

  // Generate a vector of the right type with random values
  auto rand_obj = [&]() {
    vector_soa_t tmp{};

    for (std::size_t i = 0u; i < detray::traits::size<vector_soa_t>; ++i) {
      tmp[i] = rand_simd();
    }

    return tmp;
  };

  collection.resize(collection.capacity());
  std::ranges::generate(collection, rand_obj);
}

/// Fill a lane bundle based transform3 with random values
template <detray::concepts::transform3D transform3_t>
inline void fill_random_trf(std::vector<transform3_t> &collection) {
  using vector_t = typename transform3_t::vector3;
  using simd_vector_t = typename transform3_t::scalar_type;
  using value_t = typename simd_vector_t::value_type;

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<value_t> dist(0.f, 1.f);

  auto rand_simd = [&]() {
    simd_vector_t s;
    for (std::size_t i = 0u; i < simd_vector_t::size(); ++i) {
      s[i] = dist(mt);
    }
    return s;
  };

  // Generate a random, but valid affine transformation
  auto rand_obj = [&]() {
    vector_t x_axis;
    vector_t z_axis;
    vector_t t;

    x_axis[0] = rand_simd();
    x_axis[1] = rand_simd();
    x_axis[2] = rand_simd();
    x_axis = detray::vector::normalize(x_axis);

    z_axis[0] = rand_simd();
    z_axis[1] = rand_simd();
    z_axis[2] = rand_simd();

    t[0] = rand_simd();
    t[1] = rand_simd();
    t[2] = rand_simd();
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

/// Fill a lane bundle based matrix with random values
template <detray::concepts::matrix matrix_t>
inline void fill_random_matrix(std::vector<matrix_t> &collection) {
  using simd_vector_t = typename matrix_t::scalar_type;
  using value_t = typename simd_vector_t::value_type;

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<value_t> dist(0.f, 1.f);

  auto rand_simd = [&]() {
    simd_vector_t s;
    for (std::size_t i = 0u; i < simd_vector_t::size(); ++i) {
      s[i] = dist(mt);
    }
    return s;
  };

  // Generate a random, but valid affine transformation
  auto rand_obj = [&]() {
    matrix_t m;

    for (std::size_t j = 0u; j < matrix_t::columns(); ++j) {
      typename matrix_t::vector_type v;

      for (std::size_t i = 0u; i < matrix_t::rows(); ++i) {
        v[i] = rand_simd();
      }

      m[j] = v;
    }

    return m;
  };

  collection.resize(collection.capacity());
  std::ranges::generate(collection, rand_obj);
}

}  // namespace detray::algebra
