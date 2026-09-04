// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Direction.hpp"
#include "ActsPlugins/Json/ActsJson.hpp"

#include <nlohmann/json.hpp>

namespace Acts {

/// @addtogroup json_plugin
/// @{

/// Convert Direction to JSON
/// @param j Destination JSON object
/// @param direction Source Direction to convert
void to_json(nlohmann::json& j, const Direction& direction);

/// Convert JSON to Direction
/// @param j Source JSON object
/// @param direction Destination Direction to populate
void from_json(const nlohmann::json& j, Direction& direction);

/// @}
}  // namespace Acts
