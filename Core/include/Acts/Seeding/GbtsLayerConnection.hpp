// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/GbtsLayerDescription.hpp"

#include <cstdint>
#include <vector>

namespace Acts::Experimental {

/// A pair of layers the seeder may build a graph edge between. Edges are made
/// outside-in, so a hit on @c src is the outer end of the doublet and a hit on
/// @c dst the inner one.
struct GbtsLayerConnection {
  /// Outer layer id.
  GbtsExperimentLayerId src{};
  /// Inner layer id.
  GbtsExperimentLayerId dst{};
};

/// The layers the GBTS seeding may connect
struct GbtsConnectionsConfig final {
  /// Width of the eta bins the layers are split into
  float etaBinWidth{};
  /// Pairs of connected layers
  std::vector<GbtsLayerConnection> connections;
};

}  // namespace Acts::Experimental
