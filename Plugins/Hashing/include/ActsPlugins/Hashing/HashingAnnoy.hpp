// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsPlugins/Hashing/AnnoyForwardDeclarations.hpp"

#include <map>
#include <set>

namespace ActsPlugins {
/// @addtogroup hashing_plugin
/// @{

/// Annoy-based hashing implementation for spacepoint bucketing
template <typename external_spacepoint_t, typename SpacePointContainer>
class HashingAnnoy {
 public:
  /// Compute space point buckets using Annoy model
  /// @param annoyModel Trained Annoy model
  /// @param spacePoints Space points to bucket
  /// @param bucketSize Size of each bucket
  /// @param zBins Number of z bins
  /// @param phiBins Number of phi bins
  /// @param layerRMin Minimum radial extent
  /// @param layerRMax Maximum radial extent
  /// @param layerZMin Minimum z extent
  /// @param layerZMax Maximum z extent
  void computeSpacePointsBuckets(
      const AnnoyModel* annoyModel, const SpacePointContainer& spacePoints,
      const unsigned int bucketSize, const unsigned int zBins,
      const unsigned int phiBins, const double layerRMin,
      const double layerRMax, const double layerZMin, const double layerZMax);
  /// Map of bucket indices to space points
  std::map<unsigned int, std::set<external_spacepoint_t>> m_bucketsSPMap;
};

/// @}
}  // namespace ActsPlugins

#include "ActsPlugins/Hashing/HashingAnnoy.ipp"
