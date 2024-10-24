// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/EventData/Cluster.hpp"

#include <cstdint>
#include <cstdlib>
#include <map>
#include <vector>

#include "traccc/edm/cell.hpp"

namespace ActsExamples::Traccc::Common::Conversion {

/// @brief Converts a "geometry ID -> generic cell collection type" map to a "geometry ID -> traccc cell collection" map.
/// @note The function sets the module link of the cells in the output to 0.
/// @return Map from geometry ID to its cell data (as a vector of traccc cell data)
std::map<Acts::GeometryIdentifier, std::vector<traccc::cell>> tracccCellsMap(
    const std::map<Acts::GeometryIdentifier, std::vector<Cluster::Cell>>& map);

}  // namespace ActsExamples::Traccc::Common::Conversion
