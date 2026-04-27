// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/MultiLayerNavigationJsonConverter.hpp"

#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "ActsPlugins/Json/GridJsonConverter.hpp"
#include "ActsPlugins/Json/UtilitiesJsonConverter.hpp"

nlohmann::json
Acts::MultiLayerNavigationJsonConverter::encodeMultiLayerNavigationPolicy(
    const Experimental::MultiLayerNavigationPolicy& policy,
    const TrackingGeometryJsonConverter& /*converter*/) {
  nlohmann::json jPolicy;
  jPolicy["kind"] = "MultiLayerNavigation";

  const auto& grid = policy.indexedGrid();
  nlohmann::json jAxes;
  for (const auto* axis : grid.grid.axes()) {
    jAxes.push_back(Acts::AxisJsonConverter::toJson(*axis));
  }
  jPolicy["axes"] = jAxes;

  const auto& casts = grid.casts;
  jPolicy["casts"] =
      std::vector<Acts::AxisDirection>(casts.begin(), casts.end());
  jPolicy["binExpansion"] = policy.binExpansion();

  return jPolicy;
}

std::unique_ptr<Acts::INavigationPolicy>
Acts::MultiLayerNavigationJsonConverter::decodeMultiLayerNavigationPolicy(
    const nlohmann::json& encoded, const GeometryContext& gctx,
    const TrackingGeometryJsonConverter& /*converter*/,
    const TrackingVolume& volume, const Logger& logger) {
  const auto& jAxes = encoded.at("axes");
  std::array<double, 2> range0 = jAxes.at(0).at("range");
  std::size_t bins0 = jAxes.at(0).at("bins");
  std::array<double, 2> range1 = jAxes.at(1).at("range");
  std::size_t bins1 = jAxes.at(1).at("bins");

  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axis0(range0[0],
                                                             range0[1], bins0);
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axis1(range1[0],
                                                             range1[1], bins1);
  Experimental::MultiLayerNavigationPolicy::GridType grid(std::move(axis0),
                                                          std::move(axis1));

  std::vector<AxisDirection> castsVec =
      encoded.at("casts").get<std::vector<AxisDirection>>();
  std::array<AxisDirection, 2> casts = {castsVec.at(0), castsVec.at(1)};

  Experimental::MultiLayerNavigationPolicy::IndexedUpdatorType indexedGrid(
      std::move(grid), casts);

  Experimental::MultiLayerNavigationPolicy::Config config;
  config.binExpansion =
      encoded.at("binExpansion").get<std::vector<std::size_t>>();

  return std::make_unique<Experimental::MultiLayerNavigationPolicy>(
      gctx, volume, logger, config, std::move(indexedGrid));
}
