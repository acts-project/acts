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
#include "ActsExamples/EventData/IndexSourceLink.hpp"

#include <array>
#include <cstddef>
#include <optional>
#include <span>
#include <vector>

namespace ActsExamples {

namespace {

/// The first space point of each layer, in the given order, at most @p limit
/// of them.
std::vector<SpacePointIndex> onePerLayer(
    const SpacePointContainer& spacePoints,
    std::span<const SpacePointIndex> candidates, std::size_t limit) {
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
    layers.push_back(layer);
    perLayer.push_back(index);
  }
  return perLayer;
}

}  // namespace

std::optional<std::array<SpacePointIndex, 3>> selectSeedSpacePoints(
    const SpacePointContainer& spacePoints,
    std::span<const SpacePointIndex> candidates,
    SeedSpacePointSelection selection) {
  if (candidates.size() < 3) {
    return std::nullopt;
  }

  switch (selection) {
    case SeedSpacePointSelection::FirstThree:
      return std::array{candidates[0], candidates[1], candidates[2]};
    case SeedSpacePointSelection::InnermostTriplet: {
      const std::vector<SpacePointIndex> perLayer =
          onePerLayer(spacePoints, candidates, 3);
      if (perLayer.size() < 3) {
        return std::nullopt;
      }
      return std::array{perLayer[0], perLayer[1], perLayer[2]};
    }
    case SeedSpacePointSelection::SpreadTriplet: {
      const std::vector<SpacePointIndex> perLayer =
          onePerLayer(spacePoints, candidates, candidates.size());
      if (perLayer.size() < 3) {
        return std::nullopt;
      }
      return std::array{perLayer.front(), perLayer[perLayer.size() / 2],
                        perLayer.back()};
    }
  }

  return std::nullopt;
}

}  // namespace ActsExamples
