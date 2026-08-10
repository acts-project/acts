// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/AxisFactoryJsonConverter.hpp"

#include "Acts/Utilities/AxisDefinitions.hpp"
#include "ActsPlugins/Json/GridJsonConverter.hpp"
#include "ActsPlugins/Json/UtilitiesJsonConverter.hpp"

#include <array>
#include <stdexcept>
#include <vector>

nlohmann::json Acts::AxisFactoryJsonConverter::toJson(
    const AxisFactory& axisFactory) {
  nlohmann::json j;

  std::optional<double> min;
  std::optional<double> max;
  if (axisFactory.isEquidistant()) {
    const auto& params = axisFactory.asEquidistant();
    j["type"] = "equidistant";
    j["bins"] = params.nBins;
    min = params.min;
    max = params.max;
  } else if (axisFactory.isDeferredVariable()) {
    j["type"] = "deferred-variable";
    j["normalized_boundaries"] =
        axisFactory.asDeferredVariable().normalizedEdges;
  } else {
    const auto& params = axisFactory.asVariable();
    j["type"] = "variable";
    j["boundaries"] = params.edges;
  }

  // Only the properties the description fixes are written out
  if (min.has_value() && max.has_value()) {
    j["range"] = std::array<double, 2u>({*min, *max});
  }
  if (axisFactory.boundaryType().has_value()) {
    j["boundary_type"] = *axisFactory.boundaryType();
  }
  if (axisFactory.direction().has_value()) {
    j["direction"] = *axisFactory.direction();
  }
  return j;
}

Acts::AxisFactory Acts::AxisFactoryJsonConverter::fromJson(
    const nlohmann::json& j) {
  std::optional<AxisBoundaryType> boundaryType;
  if (j.contains("boundary_type")) {
    boundaryType = j.at("boundary_type").get<AxisBoundaryType>();
  }
  std::optional<AxisDirection> direction;
  if (j.contains("direction")) {
    direction = j.at("direction").get<AxisDirection>();
  }

  std::string type = j.at("type").get<std::string>();
  if (type == "equidistant") {
    std::optional<double> min;
    std::optional<double> max;
    if (j.contains("range")) {
      std::array<double, 2u> range = j.at("range");
      min = range.at(0);
      max = range.at(1);
    }
    return AxisFactory::Equidistant(j.at("bins").get<std::size_t>(), min, max,
                                    boundaryType, direction);
  }
  if (type == "variable") {
    return AxisFactory::Variable(j.at("boundaries").get<std::vector<double>>(),
                                 boundaryType, direction);
  }
  if (type == "deferred-variable") {
    return AxisFactory::DeferredVariable(
        j.at("normalized_boundaries").get<std::vector<double>>(), boundaryType,
        direction);
  }
  throw std::invalid_argument(
      "AxisFactoryJsonConverter: unknown axis description type '" + type + "'");
}

nlohmann::json Acts::MultiAxisFactoryJsonConverter::toJson(
    const MultiAxisFactory& multiAxisFactory) {
  nlohmann::json j = nlohmann::json::array();
  for (const AxisFactory& axisFactory : multiAxisFactory.axisFactories()) {
    j.push_back(AxisFactoryJsonConverter::toJson(axisFactory));
  }
  return j;
}

Acts::MultiAxisFactory Acts::MultiAxisFactoryJsonConverter::fromJson(
    const nlohmann::json& j) {
  std::vector<AxisFactory> axisFactories;
  axisFactories.reserve(j.size());
  for (const auto& jAxis : j) {
    axisFactories.push_back(AxisFactoryJsonConverter::fromJson(jAxis));
  }
  return MultiAxisFactory(std::move(axisFactories));
}
