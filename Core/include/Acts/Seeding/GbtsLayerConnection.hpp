// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/GbtsLayerDescription.hpp"

#include <span>

namespace Acts::Experimental {

/// A pair of layers the seeder may build a graph edge between. Edges are made
/// outside-in, so a hit on @c src is the outer end of the doublet and a hit on
/// @c dst the inner one.
///
/// A triplet is cut at the larger of the tolerances its inner and outer pair
/// carry for that role.
struct GbtsLayerConnection {
  /// Outer layer id.
  GbtsExperimentLayerId src{};
  /// Inner layer id.
  GbtsExperimentLayerId dst{};
  /// Largest |tau ratio| if this is the inner pair of a triplet.
  float tauRatioCutInner = 0.007f;
  /// Largest |tau ratio| if this is the outer pair of a triplet.
  float tauRatioCutOuter = 0.007f;
};

/// The tau ratio tolerances ATLAS tuned GBTS with.
struct GbtsAdaptiveTauRatioCuts {
  /// Tolerance of a triplet nothing widens.
  float tauRatioCut = 0.007f;
  /// Widening for a pixel barrel layer skipped, or for the pixel barrel left.
  float tauRatioCorr = 0.006f;
  /// Widening, on top, for a strip layer in the triplet.
  float tauRatioCorrStrip = 0.03f;
};

/// Fill the tau ratio tolerances of the layer pairs by the ATLAS rule.
///
/// Same as the per-triplet rule, except that a triplet with its inner node off
/// the pixel barrel and its middle node on it can be widened by its outer pair.
///
/// @param connections The layer pairs, whose tolerances are overwritten
/// @param layers The layers the geometry is built from
/// @param cuts The tolerances to assign
void assignAdaptiveTauRatioCuts(std::span<GbtsLayerConnection> connections,
                                std::span<const GbtsLayerDescription> layers,
                                const GbtsAdaptiveTauRatioCuts& cuts = {});

}  // namespace Acts::Experimental
