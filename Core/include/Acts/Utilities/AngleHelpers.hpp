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

template <>
struct EtaThetaConversionTraits<float> {
  static constexpr float minTheta = 1e-6f;
  static constexpr float maxTheta = std::numbers::pi_v<float> - 1e-6f;

  static constexpr float minEta = -80.0f;
  static constexpr float maxEta = 80.0f;
};

template <>
struct EtaThetaConversionTraits<double> {
  static constexpr double minTheta = 1e-12;
  static constexpr double maxTheta = std::numbers::pi - 1e-12;

  static constexpr double minEta = -700.0;
  static constexpr double maxEta = 700.0;
};

template <typename Scalar>
struct ClampedEtaThetaConversionTraits {
  // computed from maxEta
  static constexpr Scalar minTheta =
      static_cast<Scalar>(9.0799859462585550026e-05L);
  // computed from minEta
  static constexpr Scalar maxTheta =
      static_cast<Scalar>(3.1415018537303306529L);

  static constexpr Scalar minEta = static_cast<Scalar>(-10.0);
  static constexpr Scalar maxEta = static_cast<Scalar>(10.0);

  static_assert(minTheta > static_cast<Scalar>(0),
                "minTheta must be bigger than 0");
  static_assert(maxTheta < std::numbers::pi_v<Scalar>,
                "maxTheta must be smaller than pi");
  static_assert(minEta > -std::numeric_limits<Scalar>::infinity(),
                "minEta must be bigger than -inf");
  static_assert(maxEta < std::numeric_limits<Scalar>::infinity(),
                "maxEta must be smaller than inf");
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

/// Calculate the pseudorapidity from the polar angle theta.
///
/// This function aims to be FPE safe and returns clamed eta values for theta
/// values outside the defined range.
///
/// @param theta is the polar angle in radian towards the z-axis.
///
/// @return the pseudorapidity towards the z-axis.
template <typename Scalar>
Scalar etaFromThetaClamped(Scalar theta) {
  if (theta <= ClampedEtaThetaConversionTraits<Scalar>::minTheta) {
    return ClampedEtaThetaConversionTraits<Scalar>::maxEta;
  }
  if (theta >= ClampedEtaThetaConversionTraits<Scalar>::maxTheta) {
    return ClampedEtaThetaConversionTraits<Scalar>::minEta;
  }

  return -std::log(std::tan(theta / 2));
}

/// Calculate the polar angle theta from the pseudorapidity.
///
/// This function aims to be FPE safe and returns clamed theta values for eta
/// values outside the defined range.
///
/// @param eta is the pseudorapidity towards the z-axis.
///
/// @return the polar angle in radian towards the z-axis.
template <typename Scalar>
Scalar thetaFromEtaClamped(Scalar eta) {
  if (eta <= ClampedEtaThetaConversionTraits<Scalar>::minEta) {
    return ClampedEtaThetaConversionTraits<Scalar>::maxTheta;
  }
  if (eta >= ClampedEtaThetaConversionTraits<Scalar>::maxEta) {
    return ClampedEtaThetaConversionTraits<Scalar>::minTheta;
  }

  return 2 * std::atan(std::exp(-eta));
}

}  // namespace Acts::AngleHelpers
