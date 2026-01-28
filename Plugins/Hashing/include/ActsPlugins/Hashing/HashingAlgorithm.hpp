// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/EventData/Types.hpp"
#include "ActsPlugins/Hashing/HashingModel.hpp"

#include <cstdint>

namespace ActsPlugins {
/// @addtogroup hashing_plugin
/// @{

class HashingAlgorithm {
 public:
  struct Config {
    /// Size of the buckets = number of spacepoints in the bucket
    std::uint32_t bucketSize = 10;
    /// Number of zBins
    std::uint32_t zBins = 0;
    /// Number of phiBins
    std::uint32_t phiBins = 50;

    /// Layer selection
    double layerRMin = 25;
    double layerRMax = 40;
    double layerZMin = -550;
    double layerZMax = 550;
  };

  explicit HashingAlgorithm(const Config& cfg);

  std::vector<std::vector<Acts::SpacePointIndex2>> execute(
      const AnnoyModel& annoyModel,
      const Acts::SpacePointContainer2& spacePoints) const;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
};

/// @}
}  // namespace ActsPlugins
