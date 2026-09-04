// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <numbers>

namespace Acts::AngleHelpers {

template <typename Scalar>
struct EtaThetaConversionTraits {};

/// Conversion limits for `float` eta/theta calculations.
///
/// Defines safe min/max theta and eta values for `float` precision.
template <>
struct EtaThetaConversionTraits<float> {
  /// Minimum theta value for float precision
  static constexpr float minTheta = 1e-6f;
  /// Maximum theta value for float precision
  static constexpr float maxTheta = std::numbers::pi_v<float> - minTheta;

  /// Maximum eta value for float precision
  static constexpr float maxEta = 80.0f;
  /// Minimum eta value for float precision
  static constexpr float minEta = -maxEta;
};

/// Conversion limits for `double` eta/theta calculations.
///
/// Defines safe min/max theta and eta values for `double` precision.
template <>
struct EtaThetaConversionTraits<double> {
  /// Minimum theta value for double precision
  static constexpr double minTheta = 1e-12;
  /// Maximum theta value for double precision
  static constexpr double maxTheta = std::numbers::pi - minTheta;

  /// Maximum eta value for double precision
  static constexpr double maxEta = 700.0;
  /// Minimum eta value for double precision
  static constexpr double minEta = -maxEta;
};

/// Calculate the pseudorapidity from the polar angle theta.
///
/// This function aims to be FPE safe and returns infinity for theta values
/// outside the floating point precision range.
///
/// @param theta is the polar angle in radian towards the z-axis.
///
/// @return the pseudorapidity towards the z-axis.
template <typename Scalar>
Scalar etaFromTheta(Scalar theta) {
  if (theta <= EtaThetaConversionTraits<Scalar>::minTheta) {
    return std::numeric_limits<Scalar>::infinity();
  }
  if (theta >= EtaThetaConversionTraits<Scalar>::maxTheta) {
    return -std::numeric_limits<Scalar>::infinity();
  }

  return -std::log(std::tan(theta / 2));
}

/// Calculate the polar angle theta from the pseudorapidity.
///
/// This function aims to be FPE safe and returns 0/pi for eta values outside
/// the floating point precision range.
///
/// @param eta is the pseudorapidity towards the z-axis.
///
/// @return the polar angle in radian towards the z-axis.
template <typename Scalar>
Scalar thetaFromEta(Scalar eta) {
  if (eta <= EtaThetaConversionTraits<Scalar>::minEta) {
    return std::numbers::pi_v<Scalar>;
  }
  if (eta >= EtaThetaConversionTraits<Scalar>::maxEta) {
    return 0;
  }

  return 2 * std::atan(std::exp(-eta));
}

}  // namespace Acts::AngleHelpers
