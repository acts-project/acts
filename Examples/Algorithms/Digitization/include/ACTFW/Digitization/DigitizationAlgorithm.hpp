// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <string>
#include <unordered_map>

#include "ACTFW/Framework/BareAlgorithm.hpp"
#include "ACTFW/Framework/RandomNumbers.hpp"
#include "Acts/Geometry/GeometryID.hpp"

namespace Acts {
class DigitizationModule;
class IdentifiedDetectorElement;
class PlanarModuleStepper;
class Surface;
class TrackingGeometry;
}  // namespace Acts

namespace FW {

/// Create planar clusters from simulation hits.
class DigitizationAlgorithm final : public BareAlgorithm {
 public:
  struct Config {
    /// Input collection of simulated hits.
    std::string inputSimulatedHits;
    /// Output collection of clusters.
    std::string outputClusters;
    /// Tracking geometry required to access global-to-local transforms.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    /// Module stepper for geometric clustering.
    std::shared_ptr<const Acts::PlanarModuleStepper> planarModuleStepper;
    /// Random numbers tool.
    std::shared_ptr<const RandomNumbers> randomNumbers;
  };

  /// Construct the digitization algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  DigitizationAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// Build clusters from input simulation hits.
  ///
  /// @param txt is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final override;

 private:
  struct Digitizable {
    const Acts::Surface* surface = nullptr;
    const Acts::IdentifiedDetectorElement* detectorElement = nullptr;
    const Acts::DigitizationModule* digitizer = nullptr;
  };

  Config m_cfg;
  /// Lookup container for all digitizable surfaces
  std::unordered_map<Acts::GeometryID, Digitizable> m_digitizables;
};

}  // namespace FW
