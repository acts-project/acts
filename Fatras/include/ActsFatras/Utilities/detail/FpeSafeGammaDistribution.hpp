// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <limits>
#include <random>

namespace ActsFatras::detail {

/// Draw random numbers from a gamma distribution.
/// Does not generate FPE underflow
///
/// Implements the same interface as the standard library distributions.
class FpeSafeGammaDistribution {
 public:
  /// Parameter struct that contains all distribution parameters.
  struct param_type {
    /// Parameters must link back to the host distribution.
    using distribution_type = FpeSafeGammaDistribution;

    /// Location parameter.
    ///
    /// Shape parameter
    double alpha = 0.5;
    /// Scale parameter.
    double theta = 1.0;

    /// Construct from parameters.
    param_type(double alpha_, double theta_) : alpha(alpha_), theta(theta_) {}
    // Explicitly defaulted construction and assignment
    param_type() = default;
    param_type(const param_type &) = default;
    param_type(param_type &&) = default;
    param_type &operator=(const param_type &) = default;
    param_type &operator=(param_type &&) = default;

    /// Parameters should be EqualityComparable
    friend bool operator==(const param_type &lhs, const param_type &rhs) {
      return (lhs.alpha == rhs.alpha) && (lhs.theta == rhs.theta);
    }
  };
  /// The type of the generated values.
  using result_type = double;

  /// Construct directly from the distribution parameters.
  FpeSafeGammaDistribution(double alpha, double theta) : m_cfg(alpha, theta) {}
  /// Construct from a parameter object.
  explicit FpeSafeGammaDistribution(const param_type &cfg) : m_cfg(cfg) {}
  // Explicitlely defaulted construction and assignment
  FpeSafeGammaDistribution() = default;
  FpeSafeGammaDistribution(const FpeSafeGammaDistribution &) = default;
  FpeSafeGammaDistribution(FpeSafeGammaDistribution &&) = default;
  FpeSafeGammaDistribution &operator=(const FpeSafeGammaDistribution &) =
      default;
  FpeSafeGammaDistribution &operator=(FpeSafeGammaDistribution &&) = default;

  /// Reset any possible internal state. Noop, since there is no internal state.
  void reset() {}
  /// Return the currently configured distribution parameters.
  param_type param() const { return m_cfg; }
  /// Set the distribution parameters.
  void param(const param_type &cfg) { m_cfg = cfg; }

  /// The minimum value the distribution generates.
  result_type min() const { return 0.; }
  /// The maximum value the distribution generates.
  result_type max() const { return std::numeric_limits<double>::infinity(); }

  /// Generate a random number from the configured gamma distribution.
  template <typename Generator>
  result_type operator()(Generator &generator) {
    return (*this)(generator, m_cfg);
  }
  /// Generate a random number from the given gamma distribution.
  template <typename Generator>
  result_type operator()(Generator &generator, const param_type &params) {
    if (params.alpha >= 1.) {
      std::gamma_distribution<double> gDist(params.alpha, params.theta);
      return gDist(generator);
    } else if (params.alpha <= 0.) {
      return 0.;
    }

    // This is from libstdc++; libcxx would give a different result
    std::gamma_distribution<double> gDistPlus1(params.alpha + 1., params.theta);
    // Get a random number using alpha+1
    const auto uPlus1 = gDistPlus1(generator);
    std::uniform_real_distribution<double> uDist(0., 1.);
    double u2 = uDist(generator);
    while (u2 == 0.) {
      u2 = uDist(generator);
    }
    // Check if this would cause underflow
    const double invAlpha = 1. / params.alpha;
    if (std::log(u2) * invAlpha + std::log(uPlus1) <
        std::log(std::numeric_limits<double>::min())) {
      // This would be underflow - get 0.0 right away
      return 0.;
    } else {
      return uPlus1 * std::pow(u2, invAlpha);
    }
  }

  /// Provide standard comparison operators
  friend bool operator==(const FpeSafeGammaDistribution &lhs,
                         const FpeSafeGammaDistribution &rhs) {
    return lhs.m_cfg == rhs.m_cfg;
  }

 private:
  param_type m_cfg;
};

}  // namespace ActsFatras::detail
