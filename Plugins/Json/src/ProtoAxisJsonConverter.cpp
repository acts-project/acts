// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/ProtoAxisJsonConverter.hpp"

#include "Acts/Plugins/Json/GridJsonConverter.hpp"
#include "Acts/Plugins/Json/UtilitiesJsonConverter.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"

nlohmann::json Acts::ProtoAxisJsonConverter::toJson(const Acts::ProtoAxis& pa) {
  nlohmann::json j;
  j["axis_dir"] = pa.getAxisDirection();
  j["axis"] = AxisJsonConverter::toJson(pa.getAxis());
  j["autorange"] = pa.isAutorange();
  return j;
}

Acts::ProtoAxis Acts::ProtoAxisJsonConverter::fromJson(
    const nlohmann::json& j) {
  auto axisDir = j.at("axis_dir").get<Acts::AxisDirection>();
  auto axisBoundaryType =
      j.at("axis").at("boundary_type").get<Acts::AxisBoundaryType>();
  if (auto axisType = j.at("axis").at("type").get<Acts::AxisType>();
      axisType == AxisType::Equidistant) {
    auto nbins = j.at("axis").at("bins").get<std::size_t>();
    if (nbins == 0) {
      throw std::invalid_argument("Number of bins must be positive");
    }

    if (j.at("autorange").get<bool>()) {
      return ProtoAxis(axisDir, axisBoundaryType, nbins);
    }
    auto min = j.at("axis").at("range").at(0).get<double>();
    auto max = j.at("axis").at("range").at(1).get<double>();
    if (min >= max) {
      throw std::invalid_argument("Invalid range: min must be less than max");
    }
    return ProtoAxis(axisDir, axisBoundaryType, min, max, nbins);
  }
  auto binEdges = j.at("axis").at("boundaries").get<std::vector<double>>();
  if (binEdges.size() < 2) {
    throw std::invalid_argument("At least two bin edges required");
  }
  if (!std::ranges::is_sorted(binEdges)) {
    throw std::invalid_argument("Bin edges must be sorted in ascending order");
  }
  return ProtoAxis(axisDir, axisBoundaryType, binEdges);
}
