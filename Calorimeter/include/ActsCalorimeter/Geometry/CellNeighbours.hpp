// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdint>
#include <vector>

namespace ActsCalorimeter {

/// An undirected adjacency edge supplied by the detector geometry.
/// Readout IDs keep their full detector-owned encoding. No adjacency is
/// inferred from their numeric difference, position, layer or readout system.
struct CellNeighbour {
  std::uint64_t first = 0;
  std::uint64_t second = 0;
};

using CellNeighbourContainer = std::vector<CellNeighbour>;

}  // namespace ActsCalorimeter
