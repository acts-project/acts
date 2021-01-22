// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "ActsExamples/Digitization/DigitizationConfig.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"

#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace Acts {
class Surface;
class TrackingGeometry;
}  // namespace Acts

namespace ActsExamples {

/// Algorithm that turns simulated hits into measurements by truth smearing.
class DigitzationAlgorithm final : public BareAlgorithm {
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
    /// The digitizers per GeometryIdentifier
    Acts::GeometryHierarchyMap<DigitizationConfig> digitizationConfigs;
  };

  /// Construct the smearing algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  DigitzationAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// Build measurement from simulation hits at input.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final override;

 private:
  // Support up to 4d measurements
  using Smearer =
      std::variant<ActsFatras::BoundParametersSmearer<RandomEngine, 1u>,
                   ActsFatras::BoundParametersSmearer<RandomEngine, 2u>,
                   ActsFatras::BoundParametersSmearer<RandomEngine, 3u>,
                   ActsFatras::BoundParametersSmearer<RandomEngine, 4u>>;

  // Nested smearer struct that holds geometric digitizer and smearing
  template <size_t kGeoDIM, size_t kSmearDIM = 0> 
  struct CombinedDigitizer {    
    const size_t kDIM = kGeoDIM + kSmearDIM;
    //GeometricDigitizationConfig geometric;
    //ActsFatras::BoundParametersSmearer<RandomEngine, kSmearDIM> smearing;
  };

  // Support max 4 digitization dimensions - either digital or smeared
  using Digitizer = 
    std::variant<CombinedDigitizer<1,0>, CombinedDigitizer<2,0>,
                 CombinedDigitizer<1,1>, CombinedDigitizer<1,2>,
                 CombinedDigitizer<1,3>, 
                 CombinedDigitizer<2,1>, CombinedDigitizer<2,2>,
                 CombinedDigitizer<0,1>, CombinedDigitizer<0,2>,
                 CombinedDigitizer<0,3>, CombinedDigitizer<0,4>>;

  Config m_cfg;
  //Acts::GeometryHierarchyMap<Digitizer> m_digitizers;

  /// Construct a fixed-size smearer from a configuration.
  template <size_t kSmearDIM>
  static Smearer makeSmearer(const SmearingConfig& cfg) {
    ActsFatras::BoundParametersSmearer<RandomEngine, kSmearDIM> impl;
    for (size_t i = 0; i < kSmearDIM; ++i) {
      impl.indices[i] = cfg.at(i).index;
      impl.smearFunctions[i] = cfg.at(i).smearFunction;
    }
    return impl;
  }
};

}  // namespace ActsExamples
