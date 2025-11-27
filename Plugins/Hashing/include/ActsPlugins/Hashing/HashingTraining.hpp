// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsPlugins/Hashing/AnnoyForwardDeclarations.hpp"
#include "ActsPlugins/Hashing/HashingTrainingConfig.hpp"

namespace ActsPlugins {
/// @addtogroup hashing_plugin
/// @{

template <typename SpacePointContainer>
class HashingTrainingAlgorithm {
 public:
  using Config = HashingTrainingConfig;

  explicit HashingTrainingAlgorithm(const Config &cfg);

  HashingTrainingAlgorithm() = default;
  HashingTrainingAlgorithm(
      const HashingTrainingAlgorithm<SpacePointContainer> &) = delete;
  HashingTrainingAlgorithm<SpacePointContainer> &operator=(
      const HashingTrainingAlgorithm<SpacePointContainer> &) = default;

  AnnoyModel execute(SpacePointContainer spacePoints) const;

  // / Get readonly access to the config parameters
  const Config &config() const { return m_cfg; }

 private:
  Config m_cfg;
};

/// @}
}  // namespace ActsPlugins

#include "ActsPlugins/Hashing/HashingTraining.ipp"
