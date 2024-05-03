// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/Hashing/HashingTrainingConfig.hpp"

#include <annoy/annoylib.h>
#include <annoy/kissrandom.h>

namespace Acts {

using AnnoyMetric = Annoy::AngularEuclidean;
using AnnoyModel =
    Annoy::AnnoyIndex<unsigned int, double, AnnoyMetric, Annoy::Kiss32Random,
                      Annoy::AnnoyIndexSingleThreadedBuildPolicy>;
template <typename SpacePointContainer>
class HashingTrainingAlgorithm {
 public:
  explicit HashingTrainingAlgorithm(const HashingTrainingConfig &cfg);

  /**
   * @brief Destroy the object.
   */
  ~HashingTrainingAlgorithm() = default;

  HashingTrainingAlgorithm() = default;
  HashingTrainingAlgorithm(
      const HashingTrainingAlgorithm<SpacePointContainer> &) = delete;
  HashingTrainingAlgorithm<SpacePointContainer> &operator=(
      const HashingTrainingAlgorithm<SpacePointContainer> &) = default;

  AnnoyModel execute(SpacePointContainer spacePoints) const;

  // / Get readonly access to the config parameters
  const Acts::HashingTrainingConfig &config() const { return m_cfg; }

 private:
  HashingTrainingConfig m_cfg;
};

}  // namespace Acts
#include "Acts/Seeding/Hashing/HashingTraining.ipp"
