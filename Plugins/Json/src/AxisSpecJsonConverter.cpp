// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/AxisSpecJsonConverter.hpp"

#include "Acts/Utilities/AxisDefinitions.hpp"
#include "ActsPlugins/Json/GridJsonConverter.hpp"
#include "ActsPlugins/Json/UtilitiesJsonConverter.hpp"

#include <array>
#include <stdexcept>
#include <vector>

nlohmann::json Acts::AxisSpecJsonConverter::toJson(const AxisSpec& axisSpec) {
  nlohmann::json j;

  std::optional<double> min;
  std::optional<double> max;
  if (axisSpec.isEquidistant()) {
    const auto& params = axisSpec.asEquidistant();
    j["type"] = "equidistant";
    j["bins"] = params.nBins;
    min = params.min;
    max = params.max;
  } else if (axisSpec.isDeferredVariable()) {
    j["type"] = "deferred-variable";
    j["normalized_boundaries"] = axisSpec.asDeferredVariable().normalizedEdges;
  } else {
    const auto& params = axisSpec.asVariable();
    j["type"] = "variable";
    j["boundaries"] = params.edges;
  }

  // Only the properties the spec fixes are written out
  if (min.has_value() && max.has_value()) {
    j["range"] = std::array<double, 2u>({*min, *max});
  }
  if (axisSpec.boundaryType().has_value()) {
    j["boundary_type"] = *axisSpec.boundaryType();
  }
  if (axisSpec.direction().has_value()) {
    j["direction"] = *axisSpec.direction();
  }
  return j;
}

Acts::AxisSpec Acts::AxisSpecJsonConverter::fromJson(const nlohmann::json& j) {
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
    return AxisSpec::Equidistant(j.at("bins").get<std::size_t>(), min, max,
                                 boundaryType, direction);
  }
  if (type == "variable") {
    return AxisSpec::Variable(j.at("boundaries").get<std::vector<double>>(),
                              boundaryType, direction);
  }
  if (type == "deferred-variable") {
    return AxisSpec::DeferredVariable(
        j.at("normalized_boundaries").get<std::vector<double>>(), boundaryType,
        direction);
  }
  throw std::invalid_argument(
      "AxisSpecJsonConverter: unknown axis spec type '" + type + "'");
}

nlohmann::json Acts::MultiAxisSpecJsonConverter::toJson(
    const MultiAxisSpec& multiAxisSpec) {
  nlohmann::json j = nlohmann::json::array();
  for (const AxisSpec& axisSpec : multiAxisSpec.axisSpecs()) {
    j.push_back(AxisSpecJsonConverter::toJson(axisSpec));
  }
  return j;
}

Acts::MultiAxisSpec Acts::MultiAxisSpecJsonConverter::fromJson(
    const nlohmann::json& j) {
  std::vector<AxisSpec> axisSpecs;
  axisSpecs.reserve(j.size());
  for (const auto& jAxis : j) {
    axisSpecs.push_back(AxisSpecJsonConverter::fromJson(jAxis));
  }
  return MultiAxisSpec(std::move(axisSpecs));
}
