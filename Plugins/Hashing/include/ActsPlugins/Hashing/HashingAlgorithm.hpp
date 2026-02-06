// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsPlugins/Hashing/AnnoyForwardDeclarations.hpp"
#include "ActsPlugins/Hashing/HashingAlgorithmConfig.hpp"

namespace ActsPlugins {
/// @addtogroup hashing_plugin
/// @{

/// Hashing-based algorithm for spacepoint grouping
template <typename external_spacepoint_t, typename SpacePointContainer>
class HashingAlgorithm {
 public:
  /// Configuration type for the hashing algorithm
  using Config = HashingAlgorithmConfig;

  /// Constructor
  /// @param cfg Configuration parameters
  explicit HashingAlgorithm(const Config& cfg);

  HashingAlgorithm() = default;

  /// Execute the hashing algorithm
  /// @param spacePoints Space points to process
  /// @param annoyModel Annoy model for nearest neighbor search
  /// @param outputCollection Output collection for grouped space points
  template <typename collection_t>
  void execute(SpacePointContainer& spacePoints, AnnoyModel* annoyModel,
               collection_t& outputCollection) const;

  /// Get readonly access to the config parameters
  /// @return Configuration parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
};

/// @}
}  // namespace ActsPlugins

#include "ActsPlugins/Hashing/HashingAlgorithm.ipp"
