// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/SpacePoint.hpp"

#include <array>
#include <optional>
#include <span>

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
};

/// Pick the three space points that seed a candidate track. The triplet
/// selections go by layer, not by space point.
///
/// @param spacePoints is the event space point container
/// @param candidates are the candidate space points, in track order
///        (innermost first); they are not sorted internally
/// @param selection picks which of them are used
/// @return the three selected space points, or nothing if the selection cannot
///         be made
std::optional<std::array<SpacePointIndex, 3>> selectSeedSpacePoints(
    const SpacePointContainer& spacePoints,
    std::span<const SpacePointIndex> candidates,
    SeedSpacePointSelection selection);

}  // namespace ActsExamples
