// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <limits>
#include <random>

namespace ActsFatras {

/// Draw random numbers from a Landau distribution.
///
/// Implements the same interface as the standard library distributions.
class LandauDistribution {
 public:
  /// Parameter struct that contains all distribution parameters.
  struct param_type {
    /// Parameters must link back to the host distribution.
    using distribution_type = LandauDistribution;

    /// Location parameter.
    ///
    /// @warning This is neither the mean nor the most probable value.
    double location = 0.0;
    /// Scale parameter.
    double scale = 1.0;

    /// Construct from parameters.
    /// @param location_ Location parameter value
    /// @param scale_ Scale parameter value
    param_type(double location_, double scale_)
        : location(location_), scale(scale_) {}
    // Explicitly defaulted construction and assignment
    param_type() = default;
    /// @brief Copy constructor
    param_type(const param_type &) = default;
    /// @brief Move constructor
    param_type(param_type &&) = default;
    /// @brief Copy assignment operator
    /// @return Reference to this parameter object
    param_type &operator=(const param_type &) = default;
    /// @brief Move assignment operator
    /// @return Reference to this parameter object
    param_type &operator=(param_type &&) = default;

    /// Parameters should be EqualityComparable
    friend bool operator==(const param_type &lhs, const param_type &rhs) {
      return (lhs.location == rhs.location) && (lhs.scale == rhs.scale);
    }
  };
  /// The type of the generated values.
  using result_type = double;

  /// Construct directly from the distribution parameters.
  /// @param location Location parameter of the distribution
  /// @param scale Scale parameter of the distribution
  LandauDistribution(double location, double scale) : m_cfg(location, scale) {}
  /// Construct from a parameter object.
  /// @param cfg Parameter configuration object
  explicit LandauDistribution(const param_type &cfg) : m_cfg(cfg) {}
  // Explicitlely defaulted construction and assignment
  LandauDistribution() = default;
  /// @brief Copy constructor
  LandauDistribution(const LandauDistribution &) = default;
  /// @brief Move constructor
  LandauDistribution(LandauDistribution &&) = default;
  /// @brief Copy assignment operator
  /// @return Reference to this distribution object
  LandauDistribution &operator=(const LandauDistribution &) = default;
  /// @brief Move assignment operator
  /// @return Reference to this distribution object
  LandauDistribution &operator=(LandauDistribution &&) = default;

  /// Reset any possible internal state. Noop, since there is no internal state.
  void reset() {}
  /// Return the currently configured distribution parameters.
  /// @return Current parameter configuration
  param_type param() const { return m_cfg; }
  /// Set the distribution parameters.
  /// @param cfg New parameter configuration to use
  void param(const param_type &cfg) { m_cfg = cfg; }

  /// The minimum value the distribution generates.
  /// @return Minimum possible value (negative infinity)
  result_type min() const { return -std::numeric_limits<double>::infinity(); }
  /// The maximum value the distribution generates.
  /// @return Maximum possible value (positive infinity)
  result_type max() const { return std::numeric_limits<double>::infinity(); }

  /// Generate a random number from the configured Landau distribution.
  /// @param generator Random number generator to use
  /// @return Random value from the Landau distribution
  template <typename Generator>
  result_type operator()(Generator &generator) {
    return (*this)(generator, m_cfg);
  }
  /// Generate a random number from the given Landau distribution.
  /// @param generator Random number generator to use
  /// @param params Distribution parameters to use for this generation
  /// @return Random value from the Landau distribution with given parameters
  template <typename Generator>
  result_type operator()(Generator &generator, const param_type &params) {
    const auto z = std::uniform_real_distribution<double>()(generator);
    return params.location + params.scale * quantile(z);
  }

  /// Provide standard comparison operators
  friend bool operator==(const LandauDistribution &lhs,
                         const LandauDistribution &rhs) {
    return lhs.m_cfg == rhs.m_cfg;
  }

 private:
  param_type m_cfg;

  static double quantile(double z);
};

}  // namespace ActsFatras
