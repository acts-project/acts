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
#include <cmath>
#include <cstddef>

namespace ActsExamples {

namespace {

/// The first space point of each layer, in the given order.
std::vector<SpacePointIndex> onePerLayer(
    const SpacePointContainer& spacePoints,
    std::span<const SpacePointIndex> candidates) {
  std::vector<SpacePointIndex> perLayer;
  std::vector<Acts::GeometryIdentifier> layers;
  for (const SpacePointIndex index : candidates) {
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

std::vector<SpacePointIndex> selectSeedSpacePoints(
    const SpacePointContainer& spacePoints,
    std::span<const SpacePointIndex> candidates,
    SeedSpacePointSelection selection) {
  if (candidates.size() < 3) {
    return {};
  }

  switch (selection) {
    case SeedSpacePointSelection::FirstThree:
      return {candidates.begin(), candidates.begin() + 3};
    case SeedSpacePointSelection::InnermostTriplet: {
      const std::vector<SpacePointIndex> perLayer =
          onePerLayer(spacePoints, candidates);
      if (perLayer.size() < 3) {
        return {};
      }
      return {perLayer.begin(), perLayer.begin() + 3};
    }
    case SeedSpacePointSelection::SpreadTriplet: {
      const std::vector<SpacePointIndex> perLayer =
          onePerLayer(spacePoints, candidates);
      if (perLayer.size() < 3) {
        return {};
      }
      return {perLayer.front(), perLayer[perLayer.size() / 2], perLayer.back()};
    }
    case SeedSpacePointSelection::All:
      return {candidates.begin(), candidates.end()};
  }

  return {};
}

}  // namespace ActsExamples
