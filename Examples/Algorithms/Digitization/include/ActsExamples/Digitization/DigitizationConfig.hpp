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

}  // namespace ActsExamples
