// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Hashing/HashingAlgorithm.hpp"

#include "Acts/EventData/Types.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "ActsPlugins/Hashing/HashingModel.hpp"

#include <numbers>
#include <set>

#include <annoy/annoylib.h>
#include <annoy/kissrandom.h>

namespace ActsPlugins {

namespace {

std::vector<std::vector<Acts::SpacePointIndex2>> computeSpacePointsBuckets(
    const AnnoyModel& annoyModel, const Acts::SpacePointContainer2& spacePoints,
    const std::size_t bucketSize, const std::size_t zBins,
    const std::size_t phiBins, const double layerRMin, const double layerRMax,
    const double layerZMin, const double layerZMax) {
  std::vector<std::set<Acts::SpacePointIndex2>> resultSets;

  std::size_t nBins = 0;
  if (zBins > 0) {
    nBins = zBins;
  } else if (phiBins > 0) {
    nBins = phiBins;
  } else {
    throw std::runtime_error("No bins defined");
  }
  resultSets.resize(nBins);

  const auto layerSelection = [&layerRMin, &layerRMax, &layerZMin, &layerZMax](
                                  const float r2, const float z) {
    const bool isInside =
        (r2 >= layerRMin * layerRMin && r2 <= layerRMax * layerRMax) &&
        (z >= layerZMin && z <= layerZMax);
    return isInside;
  };

  const auto getBinIndexZ = [&zBins, &layerZMin, &layerZMax](const float z) {
    const float binSize = (layerZMax - layerZMin) / zBins;
    return static_cast<int>((z - layerZMin + 0.5f * binSize) / binSize);
  };

  const auto getBinIndexPhi = [&phiBins](const float phi) {
    const float binSize = 2 * std::numbers::pi / phiBins;
    return static_cast<int>((phi + std::numbers::pi) / binSize);
  };

  const auto getBinIndex = [&zBins, &phiBins, &getBinIndexZ, &getBinIndexPhi](
                               const float z, const float phi) -> int {
    if (zBins > 0) {
      return getBinIndexZ(z);
    } else if (phiBins > 0) {
      return getBinIndexPhi(phi);
    } else {
      throw std::runtime_error("No bins defined");
    }
  };

  for (const auto spacePoint : spacePoints) {
    const float x = spacePoint.x();
    const float y = spacePoint.y();
    const float z = spacePoint.z();

    if (const float r2 = Acts::hypotSquare(x, y); !layerSelection(r2, z)) {
      continue;
    }

    const float phi = std::atan2(y, x);

    const int binIndex = getBinIndex(z, phi);
    if (binIndex < 0 || static_cast<std::uint32_t>(binIndex) >= nBins) {
      throw std::runtime_error("binIndex outside of bins covering");
    }

    /// Get the `bucketSize` closest spacePoints
    std::vector<std::uint32_t> neighborSpacePointIndices;
    annoyModel.get_nns_by_item(spacePoint.index(), bucketSize, -1,
                               &neighborSpacePointIndices, nullptr);
    resultSets[binIndex].insert(spacePoint.index());
    for (const auto& neighborSpacePointIndex : neighborSpacePointIndices) {
      resultSets[binIndex].insert(neighborSpacePointIndex);
    }
  }

  std::vector<std::vector<Acts::SpacePointIndex2>> result;
  result.reserve(resultSets.size());
  for (const auto& spSet : resultSets) {
    if (!spSet.empty()) {
      result.emplace_back(spSet.begin(), spSet.end());
    }
  }
  return result;
}

}  // namespace

HashingAlgorithm::HashingAlgorithm(const Config& cfg) : m_cfg(cfg) {
  if (m_cfg.bucketSize <= 0) {
    throw std::invalid_argument("Invalid bucket size");
  }
}

std::vector<std::vector<Acts::SpacePointIndex2>> HashingAlgorithm::execute(
    const AnnoyModel& annoyModel,
    const Acts::SpacePointContainer2& spacePoints) const {
  // Compute the buckets of spacepoints using the Annoy model
  std::vector<std::vector<Acts::SpacePointIndex2>> result =
      computeSpacePointsBuckets(
          annoyModel, spacePoints, m_cfg.bucketSize, m_cfg.zBins, m_cfg.phiBins,
          m_cfg.layerRMin, m_cfg.layerRMax, m_cfg.layerZMin, m_cfg.layerZMax);

  const std::size_t nBuckets = result.size();
  const std::size_t nSpacePoints = spacePoints.size();

  // Check if the number of buckets is greater than the number of space points
  if (nBuckets > nSpacePoints) {
    throw std::runtime_error("More buckets than the number of Space Points");
  }

  return result;
}

}  // namespace ActsPlugins
