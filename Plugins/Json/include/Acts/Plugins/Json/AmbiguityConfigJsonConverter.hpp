// This file is part of the Acts project.
//
// Copyright (C) 2023-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/AmbiguityResolution/ScoreBasedAmbiguityResolution.hpp"
#include "Acts/Plugins/Json/ActsJson.hpp"

namespace Acts {

using DetectorConfig = ScoreBasedAmbiguityResolution::DetectorConfig;

/// This function is used to convert the ScoreBasedAmbiguityResolutionCongig
/// from JSON to C++
std::pair<std::map<std::size_t, std::size_t>,
          std::vector<ScoreBasedAmbiguityResolution::DetectorConfig>>
from_json(const std::string&);

}  // namespace Acts
