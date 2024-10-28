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
struct EtaThetaConversionTraits {
  static constexpr Scalar minTheta = 0;
  static constexpr Scalar maxTheta = std::numbers::pi_v<Scalar>;

  static constexpr Scalar minEta = -std::numeric_limits<Scalar>::infinity();
  static constexpr Scalar maxEta = std::numeric_limits<Scalar>::infinity();
};

template <>
struct EtaThetaConversionTraits<float> {
  static constexpr float minTheta = 0;
  static constexpr float maxTheta = std::numbers::pi_v<float> - 1e-6f;

  static constexpr float minEta = -std::numeric_limits<float>::infinity();
  static constexpr float maxEta = std::numeric_limits<float>::infinity();
};

template <>
struct EtaThetaConversionTraits<double> {
  static constexpr double minTheta = 0;
  static constexpr double maxTheta = std::numbers::pi;

  static constexpr double minEta = -std::numeric_limits<double>::infinity();
  static constexpr double maxEta = std::numeric_limits<double>::infinity();
};

/// Calculate the pseudorapidity from the polar angle theta.
///
/// @param theta is the polar angle in radian towards the z-axis.
///
/// @return the pseudorapidity towards the z-axis.
template <typename Scalar>
Scalar etaFromTheta(Scalar theta) {
  if (theta < EtaThetaConversionTraits<Scalar>::minTheta) {
    return std::numeric_limits<Scalar>::infinity();
  }
  if (theta > EtaThetaConversionTraits<Scalar>::maxTheta) {
    return -std::numeric_limits<Scalar>::infinity();
  }

  return -std::log(std::tan(0.5 * theta));
}

/// Calculate the polar angle theta from the pseudorapidity.
///
/// @param eta is the pseudorapidity towards the z-axis.
///
/// @return the polar angle in radian towards the z-axis.
template <typename Scalar>
Scalar thetaFromEta(Scalar eta) {
  if (eta < EtaThetaConversionTraits<Scalar>::minEta) {
    return std::numbers::pi_v<Scalar>;
  }
  if (eta > EtaThetaConversionTraits<Scalar>::maxEta) {
    return 0;
  }

  return 2 * std::atan(std::exp(-eta));
}

}  // namespace Acts::AngleHelpers
