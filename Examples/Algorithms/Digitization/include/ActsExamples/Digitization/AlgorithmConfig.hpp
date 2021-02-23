// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/Digitization/DigitizationConfig.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"

#include <memory>
#include <string>

namespace ActsExamples {
namespace Digitization {

class AlgorithmConfig {
 public:
  AlgorithmConfig(const Options::Variables &vars)
      : AlgorithmConfig(vars,
                        Acts::GeometryHierarchyMap<DigitizationConfig>()){};

  AlgorithmConfig(const Options::Variables &vars,
                  Acts::GeometryHierarchyMap<DigitizationConfig> &&digiCfgs);

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
  /// Was the simple smearer requested
  const bool isSimpleSmearer;
  /// The digitizers per GeometryIdentifiers
  Acts::GeometryHierarchyMap<DigitizationConfig> digitizationConfigs;

 private:
  // Private initializer for SmearingAlgorithm
  void smearingConfig(const Options::Variables &vars);
};
}  // namespace Digitization
}  // namespace ActsExamples
