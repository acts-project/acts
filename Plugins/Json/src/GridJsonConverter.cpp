// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/GridJsonConverter.hpp"

#include "Acts/Utilities/IAxis.hpp"
#include "ActsPlugins/Json/UtilitiesJsonConverter.hpp"

nlohmann::json Acts::AxisJsonConverter::toJson(const IAxis& ia) {
  nlohmann::json jAxis;

  jAxis["boundary_type"] = ia.getBoundaryType();
  if (ia.getDirection().has_value()) {
    jAxis["direction"] = ia.getDirection().value();
  }
  // type, range, bins or boundaries
  if (ia.isEquidistant()) {
    jAxis["type"] = AxisType::Equidistant;
    jAxis["range"] = std::array<double, 2u>({ia.getMin(), ia.getMax()});
    jAxis["bins"] = ia.getNBins();
  } else {
    jAxis["type"] = AxisType::Variable;
    jAxis["boundaries"] = ia.getBinEdges();
  }
  return jAxis;
}

std::unique_ptr<Acts::IAxis> Acts::AxisJsonConverter::fromJson(
    const nlohmann::json& jAxis) {
  Acts::AxisType axisType = jAxis.at("type");
  Acts::AxisBoundaryType boundaryType = jAxis.at("boundary_type");

  std::optional<Acts::AxisDirection> direction = std::nullopt;
  if (jAxis.contains("direction")) {
    direction = jAxis.at("direction").get<Acts::AxisDirection>();
  }

  if (axisType == Acts::AxisType::Equidistant) {
    std::array<double, 2u> range = jAxis.at("range");
    return Acts::IAxis::createEquidistant(
        boundaryType, range.at(0), range.at(1), jAxis.at("bins"), direction);
  }

  return Acts::IAxis::createVariable(
      boundaryType, jAxis.at("boundaries").get<std::vector<double>>(),
      direction);
}
