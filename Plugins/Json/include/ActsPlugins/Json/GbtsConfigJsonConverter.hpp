// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/GbtsLayerConnection.hpp"
#include "Acts/Seeding/GbtsLayerConnectionTool.hpp"
#include "Acts/Seeding/GbtsLayerDescription.hpp"
#include "Acts/Seeding/detail/GbtsGraphTypes.hpp"
#include "ActsPlugins/Json/ActsJson.hpp"
#include "ActsPlugins/Json/GeometryIdentifierJsonConverter.hpp"

#include <nlohmann/json.hpp>

namespace Acts::Experimental {

/// @addtogroup json_plugin
/// @{

NLOHMANN_JSON_SERIALIZE_ENUM(GbtsLayerType, {{GbtsLayerType::Barrel, "barrel"},
                                             {GbtsLayerType::Endcap, "endcap"}})

NLOHMANN_JSON_SERIALIZE_ENUM(GbtsLayerTechnology,
                             {{GbtsLayerTechnology::Pixel, "pixel"},
                              {GbtsLayerTechnology::Strip, "strip"}})

/// Convert GbtsLayerConnection to JSON
/// @param j Destination JSON object
/// @param connection Source GbtsLayerConnection to convert
void to_json(nlohmann::json& j, const GbtsLayerConnection& connection);

/// Convert JSON to GbtsLayerConnection
/// @param j Source JSON object
/// @param connection Destination GbtsLayerConnection to populate
void from_json(const nlohmann::json& j, GbtsLayerConnection& connection);

/// Convert JSON to GbtsLayerConfig
/// @param j Source JSON object
/// @param layer Destination GbtsLayerConfig to populate
void from_json(const nlohmann::json& j, GbtsLayerConfig& layer);

/// Convert GbtsConnectionsConfig to JSON
/// @param j Destination JSON object
/// @param config Source GbtsConnectionsConfig to convert
void to_json(nlohmann::json& j, const GbtsConnectionsConfig& config);

/// Convert JSON to GbtsConnectionsConfig
/// @param j Source JSON object
/// @param config Destination GbtsConnectionsConfig to populate
void from_json(const nlohmann::json& j, GbtsConnectionsConfig& config);

/// Convert JSON to GbtsLayerConnectionTool::LayerDescription
/// @param j Source JSON object
/// @param layer Destination LayerDescription to populate
void from_json(const nlohmann::json& j,
               GbtsLayerConnectionTool::LayerDescription& layer);

namespace detail {

/// Convert JSON to GbtsTauBounds
/// @param j Source JSON object
/// @param bounds Destination GbtsTauBounds to populate
void from_json(const nlohmann::json& j, GbtsTauBounds& bounds);

}  // namespace detail

/// @}

}  // namespace Acts::Experimental
