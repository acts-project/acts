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

/// Type describing the conditions configuration of a detector module
struct conditions_data_config {
  vector2 shift{0.f, 0.f};
};

using conditions_config = Acts::GeometryHierarchyMap<conditions_data_config>;

}  // namespace traccc
