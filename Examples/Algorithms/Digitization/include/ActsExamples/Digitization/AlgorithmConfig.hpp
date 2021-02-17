// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Digitization/DigitizationConfig.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"

#include <memory>
#include <string>

namespace ActsExamples {
namespace Digitization {
  struct AlgorithmConfig {
    /// Input collection of simulated hits.
    std::string inputSimHits;
    /// Output source links collection.
    std::string outputSourceLinks;
    /// Output measurements collection.
    std::string outputMeasurements;
    /// Output cluster collection.
    std::string outputClusters;
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
    /// Was the simple smearer requested
    bool isSimpleSmearer;
  };
}  // namespace Digitization
}  // namespace ActsExamples
