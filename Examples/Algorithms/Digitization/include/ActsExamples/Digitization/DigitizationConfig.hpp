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

/// Takes as an argument the position, and a random engine
///  @return drift direction in local 3D coordinates
using DriftGenerator =
    std::function<Acts::Vector3(const Acts::Vector3&, RandomEngine&)>;
/// Takes as an argument the path length, the drift length, and a random engine
/// @return a charge to which the threshold can be applied
using ChargeGenerator = std::function<Acts::ActsScalar(
    Acts::ActsScalar, Acts::ActsScalar, RandomEngine&)>;
/// Takes as an argument the clsuter size and an random engine
/// @return a vector of uncorrelated covariance values
using VarianceGenerator =
    std::function<std::vector<Acts::ActsScalar>(size_t, size_t, RandomEngine&)>;

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
  DriftGenerator drift = [](const Acts::Vector3&,
                            RandomEngine&) -> Acts::Vector3 {
    return Acts::Vector3(0., 0., 0.);
  };
  double thickness = 0.;
  /// Charge generation
  ChargeGenerator charge = [](Acts::ActsScalar path, Acts::ActsScalar,
                              RandomEngine&) -> Acts::ActsScalar {
    return path;
  };
  double threshold = 0.;
  /// Position and Covariance generation
  bool digital = false;
  VarianceGenerator variances =
      [](size_t, size_t, RandomEngine&) -> std::vector<Acts::ActsScalar> {
    return {};
  };
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
