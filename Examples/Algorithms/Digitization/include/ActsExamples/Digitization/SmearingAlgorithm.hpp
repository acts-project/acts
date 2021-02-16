// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "ActsExamples/Digitization/SmearingConfig.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"

#include <memory>
#include <string>
#include <variant>
#include <vector>

namespace Acts {
class Surface;
class TrackingGeometry;
}  // namespace Acts

namespace ActsExamples {

/// Algorithm that turns simulated hits into measurements by truth smearing.
class SmearingAlgorithm final : public BareAlgorithm {
 public:
  struct Config {
    /// Input collection of simulated hits.
    std::string inputSimHits;
    /// Output source links collection.
    std::string outputSourceLinks;
    /// Output measurements collection.
    std::string outputMeasurements;
    /// Output collection to map measured hits to contributing particles.
    std::string outputMeasurementParticlesMap;
    /// Output collection to map measured hits to simulated hits.
    std::string outputMeasurementSimHitsMap;
    /// Tracking geometry required to access global-to-local transforms.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry = nullptr;
    /// Random numbers tool.
    std::shared_ptr<const RandomNumbers> randomNumbers = nullptr;
    /// The smearers per GeometryIdentifier
    Acts::GeometryHierarchyMap<SmearingConfig> smearers;
  };

  /// Construct the smearing algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  SmearingAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// Build measurement from simulation hits at input.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final override;

 private:
  // support up to 4d measurements
  using Smearer =
      std::variant<ActsFatras::BoundParametersSmearer<RandomEngine, 1u>,
                   ActsFatras::BoundParametersSmearer<RandomEngine, 2u>,
                   ActsFatras::BoundParametersSmearer<RandomEngine, 3u>,
                   ActsFatras::BoundParametersSmearer<RandomEngine, 4u>>;

  Config m_cfg;
  Acts::GeometryHierarchyMap<Smearer> m_smearers;

  /// Construct a fixed-size smearer from a configuration.
  template <size_t kSize>
  static Smearer makeSmearer(const SmearingConfig& cfg) {
    ActsFatras::BoundParametersSmearer<RandomEngine, kSize> impl;
    for (size_t i = 0; i < kSize; ++i) {
      impl.indices[i] = cfg.at(i).index;
      impl.smearFunctions[i] = cfg.at(i).smearFunction;
    }
    return impl;
  }
};

std::vector<std::pair<Acts::GeometryIdentifier, std::vector<Acts::BoundIndices>>>
getBoundIndices(const ActsExamples::SmearingAlgorithm::Config &cfg);

}  // namespace ActsExamples
