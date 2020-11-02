// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsFatras/Digitization/UncorrelatedHitSmearer.hpp"

#include <array>
#include <memory>
#include <string>

namespace Acts {
class Surface;
class TrackingGeometry;
}  // namespace Acts

namespace ActsExamples {

/// SmearingAlgorithm that turns simulated hits into measuremetns for Fitting
///
/// Different smearing functions can be configured for a list of supported
/// smearers (see below).
class SmearingAlgorithm final : public BareAlgorithm {
 public:
  template <Acts::BoundIndices... kParameters>
  using Smearer = std::pair<ActsFatras::BoundParametersSmearer<kParameters...>,
                            std::array<ActsFatras::SmearFunction<RandomEngine>,
                                       sizeof...(kParameters)>>;

  /// Supported smears for this example are
  /// - strip type (either in loc0 or loc1)
  /// - pixel type (in loc0 and loc1)
  /// - pixel type with direction
  /// - all above in a timed flavor
  using SupportedSmearer = std::variant<
      Smearer<Acts::eBoundLoc0>, Smearer<Acts::eBoundLoc1>,
      Smearer<Acts::eBoundLoc0, Acts::eBoundLoc1>,
      Smearer<Acts::eBoundLoc0, Acts::eBoundLoc1, Acts::eBoundPhi,
              Acts::eBoundTheta>,
      Smearer<Acts::eBoundLoc0, Acts::eBoundTime>,
      Smearer<Acts::eBoundLoc1, Acts::eBoundTime>,
      Smearer<Acts::eBoundLoc0, Acts::eBoundLoc1, Acts::eBoundTime>,
      Smearer<Acts::eBoundLoc0, Acts::eBoundLoc1, Acts::eBoundPhi,
              Acts::eBoundTheta, Acts::eBoundTime>>;

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
    Acts::GeometryHierarchyMap<SupportedSmearer> smearers;
    /// flag misconfiguration
    bool configured = false;
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
  /// The configuration struct containing the smearers
  Config m_cfg;
};

}  // namespace ActsExamples
