// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsExamples/Digitization/DigitizationConfig.hpp"
#include "ActsExamples/Digitization/Smearers.hpp"
#include "ActsExamples/Digitization/SmearingConfig.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsFatras/Digitization/UncorrelatedHitSmearer.hpp"

#include <algorithm>
#include <cstddef>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

namespace Acts {
class GeometryIdentifier;
class TrackingGeometry;
}  // namespace Acts

namespace ActsExamples {

/// Configuration struct for geometric digitization
///
/// If this is defined, then the geometric digitization
/// will create clusters with cells.
/// The BinUtility defines the segmentation and which parameters
/// are defined by this.
///
struct GeometricConfig {
  // The dimensions of the measurement
  std::vector<Acts::BoundIndices> indices = {};

  // The (multidimensional) binning definition for the segmentation of the
  // sensor
  Acts::BinUtility segmentation;

  // The thickness of the sensor
  double thickness = 0.;

  /// The charge smearer
  ActsFatras::SingleParameterSmearFunction<ActsExamples::RandomEngine>
      chargeSmearer = Digitization::Exact{};

  // The threshold below an cell activation is ignored
  double threshold = 0.;

  // Wether to assume digital readout (activation is either 0 or 1)
  bool digital = false;

  /// Charge generation (configurable via the chargeSmearer)
  Acts::ActsScalar charge(Acts::ActsScalar path, RandomEngine &rng) const {
    if (not chargeSmearer) {
      return path;
    }
    auto res = chargeSmearer(path, rng);
    if (res.ok()) {
      return std::max(0.0, res->first);
    } else {
      throw std::runtime_error(res.error().message());
    }
  }

  /// Position and Covariance generation (currently not implemented)
  /// Takes as an argument the clsuter size and an random engine
  /// @return a vector of uncorrelated covariance values
  std::vector<Acts::ActsScalar> variances(size_t /*size0*/, size_t /*size1*/,
                                          RandomEngine & /*rng*/) const {
    return {};
  };

  /// Drift generation (currently not implemented)
  /// Takes as an argument the position, and a random engine
  ///  @return drift direction in local 3D coordinates
  Acts::Vector3 drift(const Acts::Vector3 & /*position*/,
                      RandomEngine & /*rng*/) const {
    return Acts::Vector3(0., 0., 0.);
  };
};

/// Configuration struct for the Digitization algorithm
///
/// It contains:
/// - optional GeometricConfig
/// - optional SmearingConfig
struct DigiComponentsConfig {
  GeometricConfig geometricDigiConfig;
  SmearingConfig smearingDigiConfig = {};
};

class DigitizationConfig {
 public:
  DigitizationConfig(bool merge, double sigma, bool commonCorner)
      : DigitizationConfig(
            merge, sigma, commonCorner,
            Acts::GeometryHierarchyMap<DigiComponentsConfig>()){};

  DigitizationConfig(
      bool doMerge, double mergeNsigma, bool mergeCommonCorner,
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
  /// Energy deposit threshold for accepting a hit
  /// For a generic readout frontend we assume 1000 e/h pairs, in Si each
  /// e/h-pair requiers on average an energy of 3.65 eV (PDG  review 2023,
  /// Table 35.10)
  double minEnergyDeposit = 1000 * 3.65 * Acts::UnitConstants::eV;
  /// The digitizers per GeometryIdentifiers
  Acts::GeometryHierarchyMap<DigiComponentsConfig> digitizationConfigs;

  std::vector<
      std::pair<Acts::GeometryIdentifier, std::vector<Acts::BoundIndices>>>
  getBoundIndices() const;
};
}  // namespace ActsExamples
