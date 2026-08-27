// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/SpacePoint.hpp"

#include <optional>
#include <span>
#include <vector>

namespace ActsExamples {

/// Which space points of a candidate track are used to seed it, i.e. which ones
/// make up the seed or feed the track parameter estimate.
enum class SeedSpacePointSelection {
  /// The first three space points, in the order they are given.
  FirstThree,
  /// The three innermost layers.
  InnermostTriplet,
  /// The innermost, the middle and the outermost layer.
  SpreadTriplet,
  /// Every space point.
  All,
};

/// Pick the space points that seed a candidate track. The triplet selections
/// pick three, go by layer, not by space point, and skip candidates closer than
/// @p minTransverseDistance to any already picked one. The distance is
/// transverse rather than radial because that is what keeps the track parameter
/// estimation well defined.
///
/// @param spacePoints is the event space point container
/// @param candidates are the candidate space points, in track order
///        (innermost first); they are not sorted internally
/// @param selection picks which of them are used
/// @param minTransverseDistance is the minimum transverse distance between the
///        selected space points. Only the triplet selections apply it;
///        @c SeedSpacePointSelection::FirstThree and
///        @c SeedSpacePointSelection::All take the candidates as given.
/// @return the selected space points, at least three of them, or nothing if
///         the selection cannot be made
std::optional<std::vector<SpacePointIndex>> selectSeedSpacePoints(
    const SpacePointContainer& spacePoints,
    std::span<const SpacePointIndex> candidates,
    SeedSpacePointSelection selection, double minTransverseDistance);

}  // namespace ActsExamples
