// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Detray/DetrayJsonConverter.hpp"

#include "ActsPlugins/Detray/DetraySurfaceGridPolicy.hpp"

#include <memory>

#include <nlohmann/json.hpp>

namespace ActsPlugins {

namespace {

nlohmann::json encodeDetraySurfaceGridPolicy(
    const DetraySurfaceGridPolicy& policy,
    const Acts::TrackingGeometryJsonConverter& /*converter*/) {
  const DetraySurfaceGridConfig& cfg = policy.config();
  nlohmann::json j;
  j["kind"] = kDetraySurfaceGridPolicyJsonKind;
  j["vertexZHalfRange"] = cfg.vertexZHalfRange;
  j["minPt"] = cfg.minPt;
  j["bField"] = cfg.bField;
  j["nVertexZ"] = cfg.nVertexZ;
  j["nCurvature"] = cfg.nCurvature;
  j["maxSampleSpacing"] = cfg.maxSampleSpacing;
  if (cfg.numberOfBinsFactor.has_value()) {
    j["numberOfBinsFactor"] = *cfg.numberOfBinsFactor;
  }
  return j;
}

std::unique_ptr<Acts::INavigationPolicy> decodeDetraySurfaceGridPolicy(
    const nlohmann::json& j, const Acts::GeometryContext& gctx,
    const Acts::TrackingGeometryJsonConverter& /*converter*/,
    const Acts::TrackingVolume& volume, const Acts::Logger& logger) {
  DetraySurfaceGridConfig cfg;
  cfg.vertexZHalfRange = j.at("vertexZHalfRange").get<double>();
  cfg.minPt = j.at("minPt").get<double>();
  cfg.bField = j.at("bField").get<double>();
  cfg.nVertexZ = j.at("nVertexZ").get<std::size_t>();
  cfg.nCurvature = j.at("nCurvature").get<std::size_t>();
  cfg.maxSampleSpacing = j.at("maxSampleSpacing").get<double>();
  cfg.numberOfBinsFactor = std::nullopt;
  if (j.contains("numberOfBinsFactor")) {
    cfg.numberOfBinsFactor = j.at("numberOfBinsFactor").get<double>();
  }
  return std::make_unique<DetraySurfaceGridPolicy>(gctx, volume, logger, cfg);
}

}  // namespace

void addDetrayNavigationPolicies(
    Acts::TrackingGeometryJsonConverter::Config& config) {
  config.encodeNavigationPolicy.registerFunction(encodeDetraySurfaceGridPolicy);
  config.decodeNavigationPolicy.registerKind(kDetraySurfaceGridPolicyJsonKind,
                                             decodeDetraySurfaceGridPolicy);
}

}  // namespace ActsPlugins
