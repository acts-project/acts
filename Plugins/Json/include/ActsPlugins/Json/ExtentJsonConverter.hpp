// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Extent.hpp"
#include "ActsPlugins/Json/ActsJson.hpp"

#include <nlohmann/json.hpp>

namespace Acts {

/// @addtogroup json_plugin
/// @{

/// Convert Extent to JSON
/// @param j Destination JSON object
/// @param e Source Extent to convert
void to_json(nlohmann::json& j, const Extent& e);

/// Convert JSON to Extent
/// @param j Source JSON object
/// @param e Destination Extent to populate
void from_json(const nlohmann::json& j, Extent& e);

/// @}
}  // namespace Acts
