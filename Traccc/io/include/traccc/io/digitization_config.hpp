// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once
// Project includes
#include "traccc/definitions/primitives.hpp"

// Acts includes
#include <Acts/Geometry/GeometryHierarchyMap.hpp>

#include <cstdint>
#include <unordered_map>
#include <vector>

namespace traccc {

/// Type describing the digitization configuration of a detector module
struct module_digitization_config {
  std::vector<std::vector<float>> bin_edges;
  unsigned char dimensions = 2;
};

/// Type describing the digitization configuration for the whole detector
using digitization_config =
    Acts::GeometryHierarchyMap<module_digitization_config>;

}  // namespace traccc
