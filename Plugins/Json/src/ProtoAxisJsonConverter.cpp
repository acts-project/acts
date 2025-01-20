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
      j["axis"]["boundary_type"].get<Acts::AxisBoundaryType>();
  auto axisType = j["axis"]["type"].get<Acts::AxisType>();
  if (axisType == AxisType::Equidistant) {
    if (j["autorange"].get<bool>()) {
      auto nbins = j["axis"]["bins"].get<std::size_t>();
      return ProtoAxis(axisDir, axisBoundaryType, nbins);
    }
    auto min = j["axis"]["range"][0].get<double>();
    auto max = j["axis"]["range"][1].get<double>();
    auto nbins = j["axis"]["bins"].get<std::size_t>();
    return ProtoAxis(axisDir, axisBoundaryType, min, max, nbins);
  }
  auto binEdges = j["axis"]["boundaries"].get<std::vector<double>>();
  return ProtoAxis(axisDir, axisBoundaryType, binEdges);
}
