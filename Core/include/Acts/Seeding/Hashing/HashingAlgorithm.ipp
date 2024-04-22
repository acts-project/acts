// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/Hashing/HashingAnnoy.hpp"

#include <vector>

namespace Acts {
// constructor
template <typename external_spacepoint_t, typename SpacePointContainer>
HashingAlgorithm<external_spacepoint_t, SpacePointContainer>::HashingAlgorithm(
    const HashingAlgorithmConfig& cfg)
    : m_cfg(cfg) {
  if (m_cfg.bucketSize <= 0) {
    throw std::invalid_argument("Invalid bucket size");
  }
}

// function to create the buckets of spacepoints.
template <typename external_spacepoint_t, typename SpacePointContainer>
void HashingAlgorithm<external_spacepoint_t, SpacePointContainer>::execute(
    SpacePointContainer& spacePoints, AnnoyModel* annoyModel,
    GenericBackInserter<SpacePointContainer> outIt) const {
  using map_t = std::map<unsigned int, std::set<external_spacepoint_t>>;

  const std::size_t nSpacePoints = spacePoints.size();

  const unsigned int bucketSize = m_cfg.bucketSize;
  const unsigned int zBins = m_cfg.zBins;
  const unsigned int phiBins = m_cfg.phiBins;

  HashingAnnoy<external_spacepoint_t, SpacePointContainer>*
      AnnoyHashingInstance =
          new HashingAnnoy<external_spacepoint_t, SpacePointContainer>();
  AnnoyHashingInstance->ComputeSpacePointsBuckets(annoyModel, spacePoints,
                                                  bucketSize, zBins, phiBins);

  map_t bucketsSPMap = AnnoyHashingInstance->m_bucketsSPMap;
  unsigned int nBuckets = (unsigned int)bucketsSPMap.size();

  if (nBuckets > nSpacePoints) {
    throw std::runtime_error("More buckets than the number of Space Points");
  }
  for (unsigned int bucketIdx = 0; bucketIdx < nBuckets; bucketIdx++) {
    typename map_t::iterator iterator = bucketsSPMap.find(bucketIdx);
    if (iterator == bucketsSPMap.end()) {
      throw std::runtime_error("Not every bucket have been found");
    }
    std::set<external_spacepoint_t> bucketSet = iterator->second;
    SpacePointContainer bucket;
    for (external_spacepoint_t spacePoint : bucketSet) {
      bucket.push_back(spacePoint);
    }
    outIt = SpacePointContainer{bucket};
  }
}

}  // namespace Acts
