// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"

#include <memory>
#include <string>
#include <unordered_map>

namespace Acts {
class DigitizationModule;
class IdentifiedDetectorElement;
class PlanarModuleStepper;
class Surface;
class TrackingGeometry;
}  // namespace Acts

namespace ActsExamples {

/// Create planar clusters from simulation hits.
class PlanarSteppingAlgorithm final : public BareAlgorithm {
 public:
  struct Config {
    /// Input collection of simulated hits.
    std::string inputSimHits;
    /// Output collection of clusters.
    std::string outputClusters;
    /// Output source links collection.
    std::string outputSourceLinks;
    /// Output measurements collection.
    std::string outputMeasurements;
    /// Output collection to map measured hits to contributing particles.
    std::string outputMeasurementParticlesMap;
    /// Output collection to map measured hits to simulated hits.
    std::string outputMeasurementSimHitsMap;
    /// Tracking geometry required to access global-to-local transforms.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    /// Module stepper for geometric clustering.
    std::shared_ptr<const Acts::PlanarModuleStepper> planarModuleStepper;
    /// Random numbers tool.
    std::shared_ptr<const RandomNumbers> randomNumbers;
  };

  /// Construct the digitization algorithm.
  ///
  /// @param config is the algorithm configuration
  /// @param level is the logging level
  PlanarSteppingAlgorithm(Config config, Acts::Logging::Level level);

  /// Build clusters from input simulation hits.
  ///
  /// @param txt is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final override;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  struct Digitizable {
    const Acts::Surface* surface = nullptr;
    const Acts::IdentifiedDetectorElement* detectorElement = nullptr;
    const Acts::DigitizationModule* digitizer = nullptr;
  };

  Config m_cfg;
  /// Lookup container for all digitizable surfaces
  std::unordered_map<Acts::GeometryIdentifier, Digitizable> m_digitizables;
};

}  // namespace ActsExamples
