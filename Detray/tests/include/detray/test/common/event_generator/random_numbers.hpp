// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// Algebra includes
#include "algebra/utils/data_generator.hpp"

// System include(s)
#include <array>
#include <limits>
#include <random>

namespace detray::detail {

/// Wrapper for CPU random number generatrion for the @c random_track_generator
template <concepts::algebra algebra_t,
          typename distribution_t =
              std::uniform_real_distribution<dvalue<algebra_t>>,
          typename engine_t = std::mt19937_64>
struct random_numbers {
  using distribution_type = distribution_t;
  using engine_type = engine_t;
  using seed_type = typename engine_t::result_type;

  using value_type = dvalue<algebra_t>;
  using scalar_type = dscalar<algebra_t>;

  std::seed_seq m_seeds;
  engine_t m_engine;

  /// Default seed
  DETRAY_HOST
  random_numbers()
      : m_seeds{random_numbers::default_seed()}, m_engine{m_seeds} {}

  /// Different seed @param s for every instance
  DETRAY_HOST
  explicit random_numbers(seed_type s) : m_seeds{s}, m_engine{m_seeds} {}

  /// More entropy in seeds from collection @param s
  DETRAY_HOST
  explicit random_numbers(const std::vector<seed_type>& s)
      : m_seeds{s.begin(), s.end()}, m_engine{m_seeds} {}

  /// Copy constructor
  DETRAY_HOST
  random_numbers(random_numbers&& other) noexcept
      : m_engine(std::move(other.m_engine)) {}

  /// Generate random numbers in a given range
  DETRAY_HOST scalar_type
  operator()(const darray<value_type, 2> range = {
                 -std::numeric_limits<value_type>::max(),
                 std::numeric_limits<value_type>::max()}) {
    const value_type min{range[0]};
    const value_type max{range[1]};
    assert(min <= max);

    // Uniform
    if constexpr (std::is_same_v<distribution_t,
                                 std::uniform_real_distribution<value_type>>) {
      return distribution_t(min, max)(m_engine);

      // Normal
    } else if constexpr (std::is_same_v<distribution_t,
                                        std::normal_distribution<value_type>>) {
      value_type mu{min + 0.5f * (max - min)};
      return distribution_t(mu, 0.5f / 3.0f * (max - min))(m_engine);
    }
  }

  /// Explicit normal distribution around a @param mean and @param stddev
  DETRAY_HOST scalar_type normal(const value_type mean,
                                 const value_type stddev) {
    if (std::is_scalar_v<scalar_type>) {
      return (stddev == static_cast<value_type>(0))
                 ? mean
                 : std::normal_distribution<value_type>(mean, stddev)(m_engine);
    }
  }

  /// 50:50 coin toss
  DETRAY_HOST std::uint8_t coin_toss() {
    return std::uniform_int_distribution<std::uint8_t>(0u, 1u)(m_engine);
  }

  /// Get the default seed of the engine
  static constexpr seed_type default_seed() { return engine_t::default_seed; }
};

}  // namespace detray::detail
