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

/// Training algorithm for hashing models
template <typename SpacePointContainer>
class HashingTrainingAlgorithm {
 public:
  /// Configuration type
  using Config = HashingTrainingConfig;

  /// Constructor
  /// @param cfg Configuration parameters
  explicit HashingTrainingAlgorithm(const Config &cfg);

  HashingTrainingAlgorithm() = default;
  HashingTrainingAlgorithm(
      const HashingTrainingAlgorithm<SpacePointContainer> &) = delete;
  /// Assignment operator
  /// @return Reference to this object
  HashingTrainingAlgorithm<SpacePointContainer> &operator=(
      const HashingTrainingAlgorithm<SpacePointContainer> &) = default;

  /// Train and return the Annoy model
  /// @param spacePoints Space points for training
  /// @return Trained Annoy model
  AnnoyModel execute(SpacePointContainer spacePoints) const;

  /// Get readonly access to the config parameters
  /// @return Configuration parameters
  const Config &config() const { return m_cfg; }

 private:
  Config m_cfg;
};

/// @}
}  // namespace ActsPlugins

#include "ActsPlugins/Hashing/HashingTraining.ipp"
