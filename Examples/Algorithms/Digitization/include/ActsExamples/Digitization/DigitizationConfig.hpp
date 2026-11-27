// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "ActsExamples/Digitization/GeometricConfig.hpp"
#include "ActsExamples/Digitization/SmearingConfig.hpp"

#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace ActsExamples {

/// Configuration struct for the Digitization algorithm
///
/// It contains:
/// - optional GeometricConfig
/// - optional SmearingConfig
struct DigiComponentsConfig {
  GeometricConfig geometricDigiConfig;
  SmearingConfig smearingDigiConfig;
};

using DigiConfigContainer = Acts::GeometryHierarchyMap<DigiComponentsConfig>;

/// Digitization configuration keyed by sensor type
///
/// Each named type holds one full digitization configuration.
/// The geometry hierarchy map assigns a tuple name to volumes, layers or
/// surfaces so all modules of one type share a single configuration.
struct DigiSensorTypesConfig {
  /// Digitization configuration per sensor type name
  std::map<std::string, DigiComponentsConfig> types;
  /// Sensor type name per geometry identifier
  Acts::GeometryHierarchyMap<std::string> assignments;
};

/// Convert a geometry-keyed configuration into the sensor-type layout
///
/// Every rile of the input becomes its own type named "rule-<index>", so the
/// result describes exactly the same digitization as the input.
///
/// @param cfg geometry-keyed digitization configuration
/// @return the equivalent sensor-type configuration
inline DigiSensorTypesConfig toSensorTypes(const DigiConfigContainer& cfg) {
  DigiSensorTypesConfig out;
  std::vector<std::pair<Acts::GeometryIdentifier, std::string>> assignments;
  assignments.reserve(cfg.size());
  for (std::size_t i = 0; i < cfg.size(); ++i) {
    std::string name = "rule-" + std::to_string(i);
    out.types.emplace(name, cfg.valueAt(i));
    assignments.emplace_back(cfg.idAt(i), std::move(name));
  }
  out.assignments =
      Acts::GeometryHierarchyMap<std::string>(std::move(assignments));
  return out;
}

}  // namespace ActsExamples
