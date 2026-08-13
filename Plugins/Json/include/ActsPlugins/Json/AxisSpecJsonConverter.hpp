// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/AxisSpec.hpp"
#include "Acts/Utilities/MultiAxisSpec.hpp"
#include "ActsPlugins/Json/ActsJson.hpp"

#include <nlohmann/json.hpp>

/// @namespace Acts::AxisSpecJsonConverter
/// @ingroup json_plugin

namespace Acts::AxisSpecJsonConverter {

/// @addtogroup json_plugin
/// @{

/// Write the AxisSpec to a json object
///
/// @param axisSpec the axis spec to be written out
/// @return JSON object representing the axis spec
nlohmann::json toJson(const AxisSpec& axisSpec);

/// Create an AxisSpec from a json object
///
/// @param j the json object to be read from
/// @return AxisSpec created from the JSON object
AxisSpec fromJson(const nlohmann::json& j);

/// @}

}  // namespace Acts::AxisSpecJsonConverter

/// @namespace Acts::MultiAxisSpecJsonConverter
/// @ingroup json_plugin

namespace Acts::MultiAxisSpecJsonConverter {

/// @addtogroup json_plugin
/// @{

/// Write the MultiAxisSpec to a json object
///
/// @param multiAxisSpec the multi-axis spec to be written out
/// @return JSON object representing the multi-axis spec
nlohmann::json toJson(const MultiAxisSpec& multiAxisSpec);

/// Create a MultiAxisSpec from a json object
///
/// @param j the json object to be read from
/// @return MultiAxisSpec created from the JSON object
MultiAxisSpec fromJson(const nlohmann::json& j);

/// @}

}  // namespace Acts::MultiAxisSpecJsonConverter
