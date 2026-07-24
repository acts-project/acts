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
  if (axisFactory.direction().has_value()) {
    j["direction"] = axisFactory.direction().value();
  }
  if (axisFactory.isDeferred()) {
    if (axisFactory.isEquidistant()) {
      j["type"] = "deferred-equidistant";
      j["bins"] = axisFactory.asDeferredEquidistant().nBins;
    } else {
      j["type"] = "deferred-variable";
      j["normalized_edges"] = axisFactory.asDeferredVariable().normalizedEdges;
    }
    return j;
  }
  if (axisFactory.isEquidistant()) {
    const auto& params = axisFactory.asEquidistant();
    j["type"] = "equidistant";
    j["boundary_type"] = params.boundaryType;
    j["range"] = std::array<double, 2u>({params.min, params.max});
    j["bins"] = params.nBins;
  } else {
    const auto& params = axisFactory.asVariable();
    j["type"] = "variable";
    j["boundary_type"] = params.boundaryType;
    j["boundaries"] = params.edges;
  }
  return j;
}

Acts::AxisFactory Acts::AxisFactoryJsonConverter::fromJson(
    const nlohmann::json& j) {
  std::optional<AxisDirection> direction = std::nullopt;
  if (j.contains("direction")) {
    direction = j.at("direction").get<AxisDirection>();
  }

  std::string type = j.at("type").get<std::string>();
  if (type == "deferred-equidistant") {
    return AxisFactory::DeferredEquidistant(j.at("bins").get<std::size_t>(),
                                            direction);
  }
  if (type == "deferred-variable") {
    return AxisFactory::DeferredVariable(
        j.at("normalized_edges").get<std::vector<double>>(), direction);
  }
  if (type == "equidistant") {
    auto boundaryType = j.at("boundary_type").get<AxisBoundaryType>();
    std::array<double, 2u> range = j.at("range");
    return AxisFactory::Equidistant(boundaryType, range.at(0), range.at(1),
                                    j.at("bins").get<std::size_t>(), direction);
  }
  if (type == "variable") {
    auto boundaryType = j.at("boundary_type").get<AxisBoundaryType>();
    return AxisFactory::Variable(
        boundaryType, j.at("boundaries").get<std::vector<double>>(), direction);
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
