// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/Hashing/HashingAlgorithmConfig.hpp"

#include <annoy/annoylib.h>
#include <annoy/kissrandom.h>

namespace Acts {

using AnnoyMetric = Annoy::AngularEuclidean;
using AnnoyModel =
    Annoy::AnnoyIndex<unsigned int, double, AnnoyMetric, Annoy::Kiss32Random,
                      Annoy::AnnoyIndexSingleThreadedBuildPolicy>;

template <typename external_spacepoint_t, typename SpacePointContainer>
class HashingAlgorithm {
 public:
  explicit HashingAlgorithm(const HashingAlgorithmConfig& cfg);

  /**
   * @brief Destroy the object.
   */
  ~HashingAlgorithm() = default;

  HashingAlgorithm() = default;
  HashingAlgorithm(const HashingAlgorithm<external_spacepoint_t,
                                          SpacePointContainer>&) = delete;
  HashingAlgorithm<external_spacepoint_t, SpacePointContainer>& operator=(
      const HashingAlgorithm<external_spacepoint_t, SpacePointContainer>&) =
      default;

  template <typename collection_t>
  void execute(SpacePointContainer& spacePoints, AnnoyModel* annoyModel,
               collection_t& outputCollection) const;

  /// Get readonly access to the config parameters
  const Acts::HashingAlgorithmConfig& config() const { return m_cfg; }

 private:
  HashingAlgorithmConfig m_cfg;
};

}  // namespace Acts
#include "Acts/Seeding/Hashing/HashingAlgorithm.ipp"
