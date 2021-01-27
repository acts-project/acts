// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "ActsExamples/Digitization/SmearingConfig.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsFatras/Digitization/UncorrelatedHitSmearer.hpp"

#include <functional>

namespace ActsExamples {

using DriftGenerator =
    std::function<Acts::Vector3(const Acts::Vector3&, RandomEngine&)>;
using ChargeGenerator = std::function<Acts::ActsScalar(
    Acts::ActsScalar, Acts::ActsScalar, RandomEngine&)>;
using CovarianceGenerator =
    std::function<std::vector<Acts::ActsScalar>(size_t, size_t, RandomEngine&)>;

// A void drift generator - no lorentz drift.
struct VoidDriftGenerator {
  /// Charge generation call operator
  ///
  /// @param pos is the position of the call (unused here)
  /// @param rng is the random engine (unused here)
  ///
  /// @return the unchanged path length
  Acts::Vector3 operator()(const Acts::Vector3& /*pos, unused*/,
                           RandomEngine& /*rng, nused*/) const {
    return Acts::Vector3{0., 0., 0.};
  }
};

/// A void charge generator - simply returns the path length.
struct VoidChargeGenerator {
  /// Charge generation call operator
  ///
  /// @param path is the input path length in the pixel
  /// @param dlength is the drift length (unused here)
  /// @param rng is the random engine (unused here)
  ///
  /// @return the unchanged path length
  Acts::ActsScalar operator()(Acts::ActsScalar path,
                              Acts::ActsScalar /*dlength, unused*/,
                              RandomEngine& /*rng, unused*/) const {
    return path;
  }
};

/// A void covariance generator - returns empty vector
struct VoidCovarianceGenerator {
  /// Charge generation call operator
  ///
  /// @param b0 is the size of the cluster in first direction
  /// @param b1 is the size of the cluster in second direction
  /// @param rng is the random engine (unused here)
  ///
  /// @return the unchanged path length
  std::vector<Acts::ActsScalar> operator()(
      size_t /*b0, unused*/, size_t /*b1, unused*/,
      RandomEngine& /*rng, unused*/) const {
    return {};
  }
};

/// Configuration struct for geometric digitization
///
/// If this is defined, then the geometric digitization
/// will create clusters with cells.
/// The BinUtility defines the segmentation and which parameters
/// are defined by this.
///
struct GeometricDigitizationConfig {
  std::vector<Acts::BoundIndices> indices = {};
  Acts::BinUtility segmentation;
  /// Drift generation
  DriftGenerator drift = VoidDriftGenerator{};
  double thickness = 0.;
  /// Charge generation
  ChargeGenerator charge = VoidChargeGenerator{};
  double threshold = 0.;
  /// Position and Covariance generation
  bool digital = false;
  CovarianceGenerator covariance = VoidCovarianceGenerator{};
};

/// Configuration struct for the Digitization algorithm
///
/// It contains:
/// - optional GeometricConfig
/// - optional SmearingConfig
struct DigitizationConfig {
  GeometricDigitizationConfig geometricDigiConfig;
  SmearingConfig smearingDigiConfig = {};
};

}  // namespace ActsExamples
