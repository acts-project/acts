// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/GbtsConfigJsonConverter.hpp"

#include <format>
#include <string>

void Acts::Experimental::to_json(nlohmann::json& j,
                                 const GbtsLayerConnection& connection) {
  j["outer"] = connection.src;
  j["inner"] = connection.dst;
}

void Acts::Experimental::from_json(const nlohmann::json& j,
                                   GbtsLayerConnection& connection) {
  connection.src = j.at("outer").get<GbtsExperimentLayerId>();
  connection.dst = j.at("inner").get<GbtsExperimentLayerId>();
}

void Acts::Experimental::from_json(const nlohmann::json& j,
                                   GbtsLayerConfig& layer) {
  layer.id = j.at("id").get<GbtsExperimentLayerId>();
  layer.type = j.at("type").get<GbtsLayerType>();
  layer.technology = j.at("technology").get<GbtsLayerTechnology>();
  layer.surfaces = j.at("surfaces").get<std::vector<GeometryIdentifier>>();
}

void Acts::Experimental::to_json(nlohmann::json& j,
                                 const GbtsConnectionsConfig& config) {
  // the shortest decimal that reads back as the same float
  j["etaBinWidth"] = std::stod(std::format("{}", config.etaBinWidth));
  j["connections"] = config.connections;
}

void Acts::Experimental::from_json(const nlohmann::json& j,
                                   GbtsConnectionsConfig& config) {
  config.etaBinWidth = j.at("etaBinWidth").get<float>();
  config.connections =
      j.at("connections").get<std::vector<GbtsLayerConnection>>();
}

void Acts::Experimental::detail::from_json(const nlohmann::json& j,
                                           GbtsTauBounds& bounds) {
  bounds.minTau = j.at("minTau").get<float>();
  bounds.maxTau = j.at("maxTau").get<float>();
  bounds.minTauNearEdge = j.at("minTauNearEdge").get<float>();
  bounds.maxTauNearEdge = j.at("maxTauNearEdge").get<float>();
}

void Acts::Experimental::from_json(
    const nlohmann::json& j, GbtsLayerConnectionTool::LayerDescription& layer) {
  layer.minR = j.at("minR").get<float>();
  layer.maxR = j.at("maxR").get<float>();
  layer.minZ = j.at("minZ").get<float>();
  layer.maxZ = j.at("maxZ").get<float>();
  layer.gbtsId = j.at("id").get<GbtsExperimentLayerId>();
}
