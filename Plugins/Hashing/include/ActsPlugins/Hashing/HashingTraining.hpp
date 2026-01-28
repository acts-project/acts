// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SpacePointContainer2.hpp"
#include "ActsPlugins/Hashing/HashingModel.hpp"

#include <cstdint>

namespace ActsPlugins {
/// @addtogroup hashing_plugin
/// @{

class HashingTraining {
 public:
  struct Config {
    /// Random seed for Annoy
    std::uint32_t annoySeed = 123456789;

    /// Number of features to use
    std::int32_t f = 1;
  };

  explicit HashingTraining(const Config &cfg);

  AnnoyModel execute(const Acts::SpacePointContainer2 &spacePoints) const;

  // / Get readonly access to the config parameters
  const Config &config() const { return m_cfg; }

 private:
  Config m_cfg;
};

/// @}
}  // namespace ActsPlugins
