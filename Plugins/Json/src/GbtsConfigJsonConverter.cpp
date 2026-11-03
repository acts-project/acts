// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/GbtsConfigJsonConverter.hpp"

#include "ActsPlugins/Json/detail/JsonIo.hpp"

void Acts::to_json(nlohmann::json& j, const GbtsLayerConnection& connection) {
  j["outer"] = connection.src;
  j["inner"] = connection.dst;
}

void Acts::from_json(const nlohmann::json& j, GbtsLayerConnection& connection) {
  connection.src = j.at("outer").get<GbtsExperimentLayerId>();
  connection.dst = j.at("inner").get<GbtsExperimentLayerId>();
}

void Acts::to_json(nlohmann::json& j, const GbtsLayerConfig& layer) {
  j["id"] = layer.id;
  j["type"] = layer.type;
  j["technology"] = layer.technology;
  j["surfaces"] = layer.surfaces;
}

void Acts::from_json(const nlohmann::json& j, GbtsLayerConfig& layer) {
  layer.id = j.at("id").get<GbtsExperimentLayerId>();
  layer.type = j.at("type").get<GbtsLayerType>();
  layer.technology = j.at("technology").get<GbtsLayerTechnology>();
  layer.surfaces = j.at("surfaces").get<std::vector<GeometryIdentifier>>();
}

void Acts::to_json(nlohmann::json& j, const GbtsTauBounds& bounds) {
  j["minTau"] = bounds.minTau;
  j["maxTau"] = bounds.maxTau;
  j["minTauNearEdge"] = bounds.minTauNearEdge;
  j["maxTauNearEdge"] = bounds.maxTauNearEdge;
}

void Acts::from_json(const nlohmann::json& j, GbtsTauBounds& bounds) {
  bounds.minTau = j.at("minTau").get<float>();
  bounds.maxTau = j.at("maxTau").get<float>();
  bounds.minTauNearEdge = j.at("minTauNearEdge").get<float>();
  bounds.maxTauNearEdge = j.at("maxTauNearEdge").get<float>();
}

std::vector<Acts::GbtsLayerConfig> Acts::readGbtsLayers(
    const std::filesystem::path& path) {
  return Acts::detail::readJsonFile(path)
      .at("layers")
      .get<std::vector<GbtsLayerConfig>>();
}

std::vector<Acts::GbtsLayerConnection> Acts::readGbtsConnections(
    const std::filesystem::path& path) {
  return Acts::detail::readJsonFile(path)
      .at("connections")
      .get<std::vector<GbtsLayerConnection>>();
}

void Acts::writeGbtsConnections(
    const std::filesystem::path& path,
    const std::vector<GbtsLayerConnection>& connections) {
  Acts::detail::writeJsonFile(
      path, nlohmann::json{{"connections", connections}}, 4, 0);
}

Acts::GbtsTauLookupTable Acts::readGbtsTauLookupTable(
    const std::filesystem::path& path) {
  return Acts::detail::readJsonFile(path)
      .at("tauLookupTable")
      .get<GbtsTauLookupTable>();
}
