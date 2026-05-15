// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/AmbiguityResolution/ScoreBasedAmbiguityResolution.hpp"

#include <map>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

namespace Acts {

/// @addtogroup json_plugin
/// @{

/// @brief Type alias for detector-specific ambiguity resolution configuration
/// @details Configuration parameters for ambiguity resolution in a specific detector component
using DetectorConfig = ScoreBasedAmbiguityResolution::DetectorConfig;
/// @brief Type alias for a pair of detector configuration and ambiguity resolution configuration
/// @details Maps detector IDs to their configurations and stores a vector of detector-specific settings
using ConfigPair =
    std::pair<std::map<std::size_t, std::size_t>, std::vector<DetectorConfig>>;
/// Convert JSON to ConfigPair
/// @param j Source JSON object
/// @param p Destination ConfigPair to populate
void from_json(const nlohmann::json& j, ConfigPair& p);

/// @}
}  // namespace Acts
