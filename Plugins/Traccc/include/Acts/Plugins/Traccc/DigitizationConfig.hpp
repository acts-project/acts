// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Utilities/BinUtility.hpp"

namespace Acts::TracccPlugin {

/// Type describing the digitization configuration of a detector module
struct ModuleDigitizationConfig {
  Acts::BinUtility segmentation;
  char dimensions = 2;
  float variance_y = 0.f;
};

/// Type describing the digitization configuration for the whole detector
using DigitizationConfig = Acts::GeometryHierarchyMap<ModuleDigitizationConfig>;

}  // namespace Acts::TracccPlugin
