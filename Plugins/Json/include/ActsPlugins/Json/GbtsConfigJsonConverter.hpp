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
#include "Acts/Seeding/GbtsTauLookupTable.hpp"
#include "ActsPlugins/Json/ActsJson.hpp"
#include "ActsPlugins/Json/GeometryIdentifierJsonConverter.hpp"

#include <filesystem>
#include <vector>

#include <nlohmann/json.hpp>

namespace Acts::Experimental {

/// @addtogroup json_plugin
/// @{

/// @cond
NLOHMANN_JSON_SERIALIZE_ENUM(GbtsLayerType, {{GbtsLayerType::Barrel, "barrel"},
                                             {GbtsLayerType::Endcap, "endcap"}})

NLOHMANN_JSON_SERIALIZE_ENUM(GbtsLayerTechnology,
                             {{GbtsLayerTechnology::Pixel, "pixel"},
                              {GbtsLayerTechnology::Strip, "strip"}})
/// @endcond

/// Convert GbtsLayerConnection to JSON
/// @param j Destination JSON object
/// @param connection Source GbtsLayerConnection to convert
void to_json(nlohmann::json& j, const GbtsLayerConnection& connection);

/// Convert JSON to GbtsLayerConnection
/// @param j Source JSON object
/// @param connection Destination GbtsLayerConnection to populate
void from_json(const nlohmann::json& j, GbtsLayerConnection& connection);

/// Convert GbtsLayerConfig to JSON
/// @param j Destination JSON object
/// @param layer Source GbtsLayerConfig to convert
void to_json(nlohmann::json& j, const GbtsLayerConfig& layer);

/// Convert JSON to GbtsLayerConfig
/// @param j Source JSON object
/// @param layer Destination GbtsLayerConfig to populate
void from_json(const nlohmann::json& j, GbtsLayerConfig& layer);

/// Convert GbtsLayerConnectionTool::LayerDescription to JSON
/// @param j Destination JSON object
/// @param layer Source LayerDescription to convert
void to_json(nlohmann::json& j,
             const GbtsLayerConnectionTool::LayerDescription& layer);

/// Convert JSON to GbtsLayerConnectionTool::LayerDescription
/// @param j Source JSON object
/// @param layer Destination LayerDescription to populate
void from_json(const nlohmann::json& j,
               GbtsLayerConnectionTool::LayerDescription& layer);

/// Convert GbtsTauBounds to JSON
/// @param j Destination JSON object
/// @param bounds Source GbtsTauBounds to convert
void to_json(nlohmann::json& j, const GbtsTauBounds& bounds);

/// Convert JSON to GbtsTauBounds
/// @param j Source JSON object
/// @param bounds Destination GbtsTauBounds to populate
void from_json(const nlohmann::json& j, GbtsTauBounds& bounds);

/// Read the GBTS layers and their surfaces from a layer map file
/// @param path The file to read, any format the JSON plugin reads
/// @return The layers of the file's `layers` entry
std::vector<GbtsLayerConfig> readGbtsLayers(const std::filesystem::path& path);

/// Read the layer geometry the layer connection training needs from a layer
/// map file
/// @param path The file to read, any format the JSON plugin reads
/// @return The layer descriptions of the file's `layers` entry
std::vector<GbtsLayerConnectionTool::LayerDescription>
readGbtsLayerDescriptions(const std::filesystem::path& path);

/// Read the layer connections from a connection table file
/// @param path The file to read, any format the JSON plugin reads
/// @return The connections of the file's `connections` entry
std::vector<GbtsLayerConnection> readGbtsConnections(
    const std::filesystem::path& path);

/// Write layer connections to a connection table file
/// @param path The file to write, whose extension selects the format
/// @param connections The connections to write as the `connections` entry
void writeGbtsConnections(const std::filesystem::path& path,
                          const std::vector<GbtsLayerConnection>& connections);

/// Read the tau lookup table of the cluster width cuts
/// @param path The file to read, any format the JSON plugin reads
/// @return The table of the file's `tauLookupTable` entry
GbtsTauLookupTable readGbtsTauLookupTable(const std::filesystem::path& path);

/// @}

}  // namespace Acts::Experimental
