// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/EventDataTransforms.hpp"

#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/SpacePoint.hpp"

#include <algorithm>
#include <vector>

ActsExamples::ProtoTrack ActsExamples::seedToProtoTrack(
    const ConstSeedProxy& seed) {
  ProtoTrack track;
  track.reserve(seed.size());
  for (const auto& sp : seed.spacePoints()) {
    for (const auto& slink : sp.sourceLinks()) {
      const auto& islink = slink.get<IndexSourceLink>();
      track.emplace_back(islink.index());
    }
  }
  return track;
}

std::optional<ActsExamples::ConstSpacePointProxy>
ActsExamples::findSpacePointForIndex(Index index,
                                     const SpacePointContainer& spacePoints) {
  auto match = [&](const ConstSpacePointProxy& sp) {
    return std::ranges::any_of(
        sp.sourceLinks(), [&](const Acts::SourceLink& sl) {
          return sl.template get<IndexSourceLink>().index() == index;
        });
  };

  auto found = std::ranges::find_if(spacePoints, match);
  if (found == spacePoints.end()) {
    return std::nullopt;
  }
  return *found;
}

ActsExamples::SeedProxy ActsExamples::protoTrackToSeed(
    const ProtoTrack& track, const SpacePointContainer& spacePoints,
    SeedContainer& seeds) {
  auto findSpacePoint = [&](Index index) -> SpacePointIndex {
    auto found = findSpacePointForIndex(index, spacePoints);
    if (!found.has_value()) {
      throw std::runtime_error("No space point found for source-link index " +
                               std::to_string(index));
    }
    return found->index();
  };

  const auto s = track.size();
  if (s < 3) {
    throw std::runtime_error(
        "Cannot convert track with less then 3 space points to seed");
  }

  std::vector<SpacePointIndex> spacePointIndices;
  spacePointIndices.reserve(track.size());

  std::transform(track.begin(), track.end(),
                 std::back_inserter(spacePointIndices), findSpacePoint);
  std::ranges::sort(spacePointIndices, {},
                    [&](SpacePointIndex i) { return spacePoints.at(i).r(); });

  // Simply use r = m*z + t and solve for r=0 to find z vertex position...
  // Probably not the textbook way to do
  const auto& frontSp = spacePoints.at(spacePointIndices.front());
  const auto& backSp = spacePoints.at(spacePointIndices.back());
  const float m = (backSp.r() - frontSp.r()) / (backSp.z() - frontSp.z());
  const float t = frontSp.r() - m * frontSp.z();
  const float vertexZ = -t / m;

  auto seed = seeds.createSeed();
  seed.assignSpacePointIndices(spacePointIndices);
  seed.vertexZ() = vertexZ;

  return seed;
}
