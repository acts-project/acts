// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/EventData/SeedSpacePointSelection.hpp"

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <optional>
#include <span>
#include <vector>

namespace ActsExamples {

namespace {

double transverseDistance(const ConstSpacePointProxy& a,
                          const ConstSpacePointProxy& b) {
  return Acts::fastHypot(a.x() - b.x(), a.y() - b.y());
}

/// The first space point of each layer, in the given order, at most @p limit of
/// them, keeping every pick at least @p minTransverseDistance away from all the
/// previous ones.
std::vector<SpacePointIndex> onePerLayer(
    const SpacePointContainer& spacePoints,
    std::span<const SpacePointIndex> candidates, std::size_t limit,
    double minTransverseDistance) {
  std::vector<SpacePointIndex> perLayer;
  std::vector<Acts::GeometryIdentifier> layers;
  for (const SpacePointIndex index : candidates) {
    if (perLayer.size() == limit) {
      break;
    }
    const ConstSpacePointProxy sp = spacePoints.at(index);
    if (sp.sourceLinks().empty()) {
      continue;
    }
    const Acts::GeometryIdentifier layer =
        sp.sourceLinks()[0].get<IndexSourceLink>().geometryId().withSensitive(
            0);
    if (Acts::rangeContainsValue(layers, layer)) {
      continue;
    }
    // a new layer does not imply a new transverse position, e.g. a forward
    // track crosses successive endcap disks at nearly the same one, and a
    // curling track comes back to where it started
    const bool tooClose =
        std::ranges::any_of(perLayer, [&](const SpacePointIndex picked) {
          return transverseDistance(sp, spacePoints.at(picked)) <
                 minTransverseDistance;
        });
    if (tooClose) {
      continue;
    }
    layers.push_back(layer);
    perLayer.push_back(index);
  }
  return perLayer;
}

}  // namespace

}  // namespace ActsExamples

std::optional<std::array<ActsExamples::SpacePointIndex, 3>>
ActsExamples::selectSeedSpacePoints(const SpacePointContainer& spacePoints,
                                    std::span<const SpacePointIndex> candidates,
                                    SeedSpacePointSelection selection,
                                    double minTransverseDistance) {
  if (candidates.size() < 3) {
    return std::nullopt;
  }

  switch (selection) {
    case SeedSpacePointSelection::FirstThree:
      return std::array{candidates[0], candidates[1], candidates[2]};
    case SeedSpacePointSelection::InnermostTriplet: {
      const std::vector<SpacePointIndex> perLayer =
          onePerLayer(spacePoints, candidates, 3, minTransverseDistance);
      if (perLayer.size() < 3) {
        return std::nullopt;
      }
      return std::array{perLayer[0], perLayer[1], perLayer[2]};
    }
    case SeedSpacePointSelection::SpreadTriplet: {
      const std::vector<SpacePointIndex> perLayer = onePerLayer(
          spacePoints, candidates, candidates.size(), minTransverseDistance);
      if (perLayer.size() < 3) {
        return std::nullopt;
      }
      return std::array{perLayer.front(), perLayer[perLayer.size() / 2],
                        perLayer.back()};
    }
    case SeedSpacePointSelection::All:
      // no triplet to pick
      return std::nullopt;
  }

  return std::nullopt;
}
