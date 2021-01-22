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
#include "ActsFatras/Digitization/UncorrelatedHitSmearer.hpp"

namespace ActsExamples {

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
  // The drift information
  double thickness = 0.;
  double threshold = 0.;
  bool digital = true;
  Acts::Vector3 driftDirection = Acts::Vector3(0., 0., 0.);
  ParameterSmearingConfig driftSmearing;
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
