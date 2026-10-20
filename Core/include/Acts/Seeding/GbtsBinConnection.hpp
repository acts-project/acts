// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/GbtsLayerDescription.hpp"

#include <compare>
#include <cstdint>

namespace Acts::Experimental {

/// One eta bin of a GBTS layer.
struct GbtsLayerBin {
  /// Id of the layer the bin belongs to.
  GbtsExperimentLayerId layer{};
  /// Index of the bin inside that layer, not in the global bin numbering.
  std::uint32_t bin{};

  /// Order two bins, by layer id and then bin index
  /// @param lhs The first bin
  /// @param rhs The second bin
  /// @return The ordering of the two bins
  friend auto operator<=>(const GbtsLayerBin& lhs,
                          const GbtsLayerBin& rhs) = default;
};

/// A pair of eta bins the seeder may build a graph edge between. Edges are
/// made outside-in, so a hit in @c src is the outer end of the doublet and a
/// hit in @c dst the inner one.
struct GbtsBinConnection {
  /// Outer bin.
  GbtsLayerBin src{};
  /// Inner bin.
  GbtsLayerBin dst{};
};

}  // namespace Acts::Experimental
