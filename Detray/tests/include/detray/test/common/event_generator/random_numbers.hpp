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

// System include(s)
#include <array>
#include <limits>
#include <random>

namespace detray::detail {

/// Wrapper for CPU random number generatrion for the @c random_track_generator
template <concepts::scalar scalar_t,
          typename distribution_t = std::uniform_real_distribution<scalar_t>,
          typename engine_t = std::mt19937_64>
struct random_numbers {
  /// Random distribution type
  using distribution_type = distribution_t;
  /// Random engine type
  using engine_type = engine_t;
  /// Random seed type
  using seed_type = typename engine_t::result_type;

  /// Seed sequence
  std::seed_seq m_seeds;
  /// Random engine
  engine_t m_engine;

  /// Default seed
  DETRAY_HOST
  random_numbers()
      : m_seeds{random_numbers::default_seed()}, m_engine{m_seeds} {}

  /// Different seed @p s for every instance
  /// @param s random engine seed
  DETRAY_HOST
  explicit random_numbers(seed_type s) : m_seeds{s}, m_engine{m_seeds} {}

  /// More entropy in seeds from collection @p s
  /// @param s random engine seeds
  DETRAY_HOST
  explicit random_numbers(const std::vector<seed_type>& s)
      : m_seeds{s.begin(), s.end()}, m_engine{m_seeds} {}

  /// Move constructor
  /// @param other random number generator to move
  DETRAY_HOST
  random_numbers(random_numbers&& other) noexcept
      : m_engine(std::move(other.m_engine)) {}

  /// Generate random numbers in a given range
  /// @param range lower and upper distribution bounds
  /// @return the generated random value
  DETRAY_HOST auto operator()(const darray<scalar_t, 2> range = {
                                  -std::numeric_limits<scalar_t>::max(),
                                  std::numeric_limits<scalar_t>::max()}) {
    const scalar_t min{range[0]};
    const scalar_t max{range[1]};
    assert(min <= max);

    // Uniform
    if constexpr (std::is_same_v<distribution_t,
                                 std::uniform_real_distribution<scalar_t>>) {
      return distribution_t(min, max)(m_engine);

      // Normal
    } else if constexpr (std::is_same_v<distribution_t,
                                        std::normal_distribution<scalar_t>>) {
      scalar_t mu{min + 0.5f * (max - min)};
      return distribution_t(mu, 0.5f / 3.0f * (max - min))(m_engine);
    }
  }

  /// Explicit normal distribution around a @p mean and @p stddev
  /// @param mean distribution mean
  /// @param stddev distribution standard deviation
  /// @return the generated random value
  DETRAY_HOST auto normal(const scalar_t mean, const scalar_t stddev) {
    return (stddev == static_cast<scalar_t>(0))
               ? mean
               : std::normal_distribution<scalar_t>(mean, stddev)(m_engine);
  }

  /// 50:50 coin toss
  /// @return zero or one
  DETRAY_HOST std::uint8_t coin_toss() {
    return std::uniform_int_distribution<std::uint8_t>(0u, 1u)(m_engine);
  }

  /// Get the default seed of the engine
  /// @return the default seed
  static constexpr seed_type default_seed() { return engine_t::default_seed; }
};

}  // namespace detray::detail
