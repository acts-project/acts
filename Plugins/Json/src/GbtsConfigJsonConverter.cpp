// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/GbtsConfigJsonConverter.hpp"

#include "ActsPlugins/Json/detail/JsonIo.hpp"

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

void Acts::Experimental::to_json(nlohmann::json& j,
                                 const GbtsLayerConfig& layer) {
  j["id"] = layer.id;
  j["type"] = layer.type;
  j["technology"] = layer.technology;
  j["surfaces"] = layer.surfaces;
}

void Acts::Experimental::from_json(const nlohmann::json& j,
                                   GbtsLayerConfig& layer) {
  layer.id = j.at("id").get<GbtsExperimentLayerId>();
  layer.type = j.at("type").get<GbtsLayerType>();
  layer.technology = j.at("technology").get<GbtsLayerTechnology>();
  layer.surfaces = j.at("surfaces").get<std::vector<GeometryIdentifier>>();
}

void Acts::Experimental::to_json(nlohmann::json& j,
                                 const GbtsTauBounds& bounds) {
  j["minTau"] = bounds.minTau;
  j["maxTau"] = bounds.maxTau;
  j["minTauNearEdge"] = bounds.minTauNearEdge;
  j["maxTauNearEdge"] = bounds.maxTauNearEdge;
}

void Acts::Experimental::from_json(const nlohmann::json& j,
                                   GbtsTauBounds& bounds) {
  bounds.minTau = j.at("minTau").get<float>();
  bounds.maxTau = j.at("maxTau").get<float>();
  bounds.minTauNearEdge = j.at("minTauNearEdge").get<float>();
  bounds.maxTauNearEdge = j.at("maxTauNearEdge").get<float>();
}

void Acts::Experimental::to_json(
    nlohmann::json& j, const GbtsLayerConnectionTool::LayerDescription& layer) {
  j["minR"] = layer.minR;
  j["maxR"] = layer.maxR;
  j["minZ"] = layer.minZ;
  j["maxZ"] = layer.maxZ;
  j["id"] = layer.gbtsId;
}

void Acts::Experimental::from_json(
    const nlohmann::json& j, GbtsLayerConnectionTool::LayerDescription& layer) {
  layer.minR = j.at("minR").get<float>();
  layer.maxR = j.at("maxR").get<float>();
  layer.minZ = j.at("minZ").get<float>();
  layer.maxZ = j.at("maxZ").get<float>();
  layer.gbtsId = j.at("id").get<GbtsExperimentLayerId>();
}

std::vector<Acts::Experimental::GbtsLayerConfig>
Acts::Experimental::readGbtsLayers(const std::filesystem::path& path) {
  return Acts::detail::readJsonFile(path)
      .at("layers")
      .get<std::vector<GbtsLayerConfig>>();
}

std::vector<Acts::Experimental::GbtsLayerConnectionTool::LayerDescription>
Acts::Experimental::readGbtsLayerDescriptions(
    const std::filesystem::path& path) {
  return Acts::detail::readJsonFile(path)
      .at("layers")
      .get<std::vector<GbtsLayerConnectionTool::LayerDescription>>();
}

std::vector<Acts::Experimental::GbtsLayerConnection>
Acts::Experimental::readGbtsConnections(const std::filesystem::path& path) {
  return Acts::detail::readJsonFile(path)
      .at("connections")
      .get<std::vector<GbtsLayerConnection>>();
}

void Acts::Experimental::writeGbtsConnections(
    const std::filesystem::path& path,
    const std::vector<GbtsLayerConnection>& connections) {
  Acts::detail::writeJsonFile(
      path, nlohmann::json{{"connections", connections}}, 4, 0);
}

Acts::Experimental::GbtsTauLookupTable
Acts::Experimental::readGbtsTauLookupTable(const std::filesystem::path& path) {
  return Acts::detail::readJsonFile(path)
      .at("tauLookupTable")
      .get<GbtsTauLookupTable>();
}
