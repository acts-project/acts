// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <random>

#include "ActsFatras/Kernel/detail/LandauQuantile.hpp"

namespace ActsFatras {

/// The following standard random number distributions are supported:
///
using GaussDist = std::normal_distribution<double>;          ///< Normal
using UniformDist = std::uniform_real_distribution<double>;  ///< Uniform
using GammaDist = std::gamma_distribution<double>;           ///< Gamma
using PoissonDist = std::poisson_distribution<int>;          ///< Poisson
///
/// In addition, the Landau distribution is provided
///
class LandauDist {
 public:
  /// A RandomNumberDistribution should provide a parameters struct
  struct param_type {
    double mean = 0.;   ///< Mean of the Landau distribution
    double scale = 1.;  ///< Scale factor

    /// Default constructor and constructor from raw parameters
    param_type() = default;
    param_type(double mean_, double scale_);

    /// Parameters should be CopyConstructible and CopyAssignable
    param_type(const param_type &) = default;
    param_type &operator=(const param_type &) = default;

    /// Parameters should be EqualityComparable
    bool operator==(const param_type &other) const;
    bool operator!=(const param_type &other) const { return !(*this == other); }

    /// Parameters should link back to the host distribution
    using distribution_type = LandauDist;
  };

  /// There should be a default constructor, a constructor from raw parameters,
  /// and a constructor from a parameters struct
  LandauDist() = default;
  LandauDist(double mean, double scale);
  LandauDist(const param_type &cfg);

  /// A distribution should be copy-constructible and copy-assignable
  LandauDist(const LandauDist &) = default;
  LandauDist &operator=(const LandauDist &) = default;

  /// Some standard ways to control the distribution's state should be provided
  void reset() { /* There is currently no state to reset here */
  }
  param_type param() const { return m_cfg; }
  void param(const param_type &p) { m_cfg = p; }

  /// A RandomNumberDistribution should provide a result type typedef and some
  /// bounds on the values that can be emitted as output
  using result_type = double;
  result_type min() const;
  result_type max() const;

  /// Generate a random number following a Landau distribution
  template <typename Generator>
  result_type operator()(Generator &engine) {
    return (*this)(engine, m_cfg);
  }

  /// Do the same, but using custom Landau distribution parameters
  template <typename Generator>
  result_type operator()(Generator &engine, const param_type &params) {
    double x = std::generate_canonical<float, 10>(engine);
    double res = params.mean + landau_quantile(x, params.scale);
    return res;
  }

  /// Provide standard comparison operators
  bool operator==(const LandauDist &other) const;
  bool operator!=(const LandauDist &other) const { return !(*this == other); }

 private:
  param_type m_cfg;  ///< configuration struct
};

}  // namespace ActsFatras
