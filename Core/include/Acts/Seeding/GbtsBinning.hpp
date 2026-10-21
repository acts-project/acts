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

namespace Acts::Experimental {

/// The eta binning `GbtsGeometry` gave one layer. Eta bins are numbered
/// globally across the geometry and a layer's are contiguous, so the layer owns
/// `[firstBin, firstBin + numBins)`.
struct GbtsLayerBinning final {
  /// Global index of the layer's first eta bin.
  std::uint32_t firstBin{};
  /// Number of eta bins the layer was split into. At least one.
  std::uint32_t numBins{1};
  /// Lower edge of the first bin, padded outwards by a small epsilon so a hit
  /// on the layer boundary still falls in a bin.
  float minEta{};
  /// Width of one bin. Divided out of the layer's eta extent, so it is not the
  /// width the geometry was asked for.
  float etaBinWidth{};
};

/// The outer eta bins reachable from one inner eta bin. Groups come in the
/// order the graph is built, outside-in.
struct GbtsBinGroup final {
  /// The inner eta bin.
  std::uint32_t bin{};
  /// The outer bins it links to.
  std::vector<std::uint32_t> links;
};

}  // namespace Acts::Experimental
