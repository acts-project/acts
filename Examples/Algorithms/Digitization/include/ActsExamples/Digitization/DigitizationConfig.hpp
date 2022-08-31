// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Digitization/DigitizationConfig.hpp"
#include "ActsExamples/Digitization/SmearingConfig.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"
#include "ActsFatras/Digitization/UncorrelatedHitSmearer.hpp"

#include <functional>
#include <memory>
#include <string>

namespace ActsExamples {

/// Takes as an argument the position, and a random engine
///  @return drift direction in local 3D coordinates
using DriftGenerator =
    std::function<Acts::Vector3(const Acts::Vector3 &, RandomEngine &)>;
/// Takes as an argument the path length, the drift length, and a random engine
/// @return a charge to which the threshold can be applied
using ChargeGenerator = std::function<Acts::ActsScalar(
    Acts::ActsScalar, Acts::ActsScalar, RandomEngine &)>;
/// Takes as an argument the clsuter size and an random engine
/// @return a vector of uncorrelated covariance values
using VarianceGenerator = std::function<std::vector<Acts::ActsScalar>(
    size_t, size_t, RandomEngine &)>;

/// Configuration struct for geometric digitization
///
/// If this is defined, then the geometric digitization
/// will create clusters with cells.
/// The BinUtility defines the segmentation and which parameters
/// are defined by this.
///
struct GeometricConfig {
  std::vector<Acts::BoundIndices> indices = {};
  Acts::BinUtility segmentation;
  /// Drift generation
  DriftGenerator drift = [](const Acts::Vector3 &,
                            RandomEngine &) -> Acts::Vector3 {
    return Acts::Vector3(0., 0., 0.);
  };
  double thickness = 0.;
  /// Charge generation
  ChargeGenerator charge = [](Acts::ActsScalar path, Acts::ActsScalar,
                              RandomEngine &) -> Acts::ActsScalar {
    return path;
  };
  double threshold = 0.;
  /// Position and Covariance generation
  bool digital = false;
  VarianceGenerator variances =
      [](size_t, size_t, RandomEngine &) -> std::vector<Acts::ActsScalar> {
    return {};
  };

  /// Equality operator for basic parameters
  /// check if the geometry config can be reused from
  /// @param other, @return a boolean to indicate this
  bool operator==(const GeometricConfig &other) const {
    return (indices == other.indices and segmentation == other.segmentation and
            thickness == other.thickness and threshold == other.threshold and
            digital == other.digital);
  }
};

/// Configuration struct for the Digitization algorithm
///
/// It contains:
/// - optional GeometricConfig
/// - optional SmearingConfig
struct DigiComponentsConfig {
  GeometricConfig geometricDigiConfig;
  SmearingConfig smearingDigiConfig = {};

  /// Equality operator to check if a digitization configuration
  /// can be reused from @param other
  ///
  /// @return a boolean flag indicating equality
  bool operator==(const DigiComponentsConfig &other) const {
    return (geometricDigiConfig == other.geometricDigiConfig and
            smearingDigiConfig == other.smearingDigiConfig);
  }
};

class DigitizationConfig {
 public:
  DigitizationConfig(const Options::Variables &vars)
      : DigitizationConfig(
            vars, Acts::GeometryHierarchyMap<DigiComponentsConfig>()){};

  DigitizationConfig(
      const Options::Variables &vars,
      Acts::GeometryHierarchyMap<DigiComponentsConfig> &&digiCfgs);

  DigitizationConfig(
      Acts::GeometryHierarchyMap<DigiComponentsConfig> &&digiCfgs);

  /// Input collection of simulated hits.
  std::string inputSimHits = "simhits";
  /// Output source links collection.
  std::string outputSourceLinks = "sourcelinks";
  /// Output measurements collection.
  std::string outputMeasurements = "measurements";
  /// Output cluster collection.
  std::string outputClusters = "clusters";
  /// Output collection to map measured hits to contributing particles.
  std::string outputMeasurementParticlesMap = "measurement_particles_map";
  /// Output collection to map measured hits to simulated hits.
  std::string outputMeasurementSimHitsMap = "measurement_simhits_map";
  /// Tracking geometry required to access global-to-local transforms.
  std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry = nullptr;
  /// Random numbers tool.
  std::shared_ptr<const RandomNumbers> randomNumbers = nullptr;
  /// Do we merge hits or not
  bool doMerge;
  /// How close do parameters have to be to consider merged
  const double mergeNsigma;
  /// Consider clusters that share a corner as merged (8-cell connectivity)
  const bool mergeCommonCorner;
  /// The digitizers per GeometryIdentifiers
  Acts::GeometryHierarchyMap<DigiComponentsConfig> digitizationConfigs;

  std::vector<
      std::pair<Acts::GeometryIdentifier, std::vector<Acts::BoundIndices>>>
  getBoundIndices() const;
};
}  // namespace ActsExamples
