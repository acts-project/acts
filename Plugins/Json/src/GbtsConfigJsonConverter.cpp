// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/GbtsConfigJsonConverter.hpp"

#include "ActsPlugins/Json/detail/JsonIo.hpp"

void Acts::Experimental::to_json(nlohmann::json& j, const GbtsLayerBin& bin) {
  j["layerId"] = bin.layer;
  j["bin"] = bin.bin;
}

void Acts::Experimental::from_json(const nlohmann::json& j, GbtsLayerBin& bin) {
  bin.layer = j.at("layerId").get<GbtsExperimentLayerId>();
  bin.bin = j.at("bin").get<std::uint32_t>();
}

void Acts::Experimental::to_json(nlohmann::json& j,
                                 const GbtsBinConnection& connection) {
  j["outer"] = connection.src;
  j["inner"] = connection.dst;
}

void Acts::Experimental::from_json(const nlohmann::json& j,
                                   GbtsBinConnection& connection) {
  connection.src = j.at("outer").get<GbtsLayerBin>();
  connection.dst = j.at("inner").get<GbtsLayerBin>();
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

std::vector<Acts::Experimental::GbtsLayerConfig>
Acts::Experimental::readGbtsLayers(const std::filesystem::path& path) {
  return Acts::detail::readJsonFile(path)
      .at("layers")
      .get<std::vector<GbtsLayerConfig>>();
}

std::vector<Acts::Experimental::GbtsBinConnection>
Acts::Experimental::readGbtsConnections(const std::filesystem::path& path) {
  return Acts::detail::readJsonFile(path)
      .at("connections")
      .get<std::vector<GbtsBinConnection>>();
}

void Acts::Experimental::writeGbtsConnections(
    const std::filesystem::path& path,
    const std::vector<GbtsBinConnection>& connections) {
  Acts::detail::writeJsonFile(
      path, nlohmann::json{{"connections", connections}}, 4, 0);
}

Acts::Experimental::GbtsTauLookupTable
Acts::Experimental::readGbtsTauLookupTable(const std::filesystem::path& path) {
  return Acts::detail::readJsonFile(path)
      .at("tauLookupTable")
      .get<GbtsTauLookupTable>();
}
