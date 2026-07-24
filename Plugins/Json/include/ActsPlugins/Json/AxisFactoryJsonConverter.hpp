// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/AxisFactory.hpp"
#include "Acts/Utilities/MultiAxisFactory.hpp"
#include "ActsPlugins/Json/ActsJson.hpp"

#include <nlohmann/json.hpp>

/// @namespace Acts::AxisFactoryJsonConverter
/// @ingroup json_plugin

namespace Acts::AxisFactoryJsonConverter {

/// @addtogroup json_plugin
/// @{

/// Write the AxisFactory to a json object
///
/// @param axisFactory the axis description to be written out
/// @return JSON object representing the axis description
nlohmann::json toJson(const AxisFactory& axisFactory);

/// Create an AxisFactory from a json object
///
/// @param j the json object to be read from
/// @return AxisFactory created from the JSON object
AxisFactory fromJson(const nlohmann::json& j);

/// @}

}  // namespace Acts::AxisFactoryJsonConverter

/// @namespace Acts::MultiAxisFactoryJsonConverter
/// @ingroup json_plugin

namespace Acts::MultiAxisFactoryJsonConverter {

/// @addtogroup json_plugin
/// @{

/// Write the MultiAxisFactory to a json object
///
/// @param multiAxisFactory the multi-axis description to be written out
/// @return JSON object representing the multi-axis description
nlohmann::json toJson(const MultiAxisFactory& multiAxisFactory);

/// Create a MultiAxisFactory from a json object
///
/// @param j the json object to be read from
/// @return MultiAxisFactory created from the JSON object
MultiAxisFactory fromJson(const nlohmann::json& j);

/// @}

}  // namespace Acts::MultiAxisFactoryJsonConverter
