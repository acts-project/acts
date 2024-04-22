// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/Hashing/Annoylib.hpp"
#include "Acts/Seeding/Hashing/Kissrandom.hpp"

#include <map>
#include <set>

namespace Acts {

template <typename external_spacepoint_t, typename SpacePointContainer>
class HashingAnnoy {
 public:
  void ComputeSpacePointsBuckets(
      const Annoy::AnnoyIndex<
          unsigned int, double, Annoy::AngularEuclidean, Annoy::Kiss32Random,
          Annoy::AnnoyIndexSingleThreadedBuildPolicy>* annoyModel,
      const SpacePointContainer& spacePoints, const unsigned int bucketSize,
      const unsigned int zBins, const unsigned int phiBins);
  std::map<unsigned int, std::set<external_spacepoint_t>> m_bucketsSPMap;
};
}  // namespace Acts
#include "Acts/Seeding/Hashing/HashingAnnoy.ipp"
