// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/UtilitiesJsonConverter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "ActsPlugins/Json/AlgebraJsonConverter.hpp"
#include "ActsPlugins/Json/GridJsonConverter.hpp"

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

Acts::ProtoAxis protoAxisFromJsonImpl(const nlohmann::json& j) {
  const auto axisBoundaryType =
      j.at("axis").at("boundary_type").get<Acts::AxisBoundaryType>();
  const auto axisType = j.at("axis").at("type").get<Acts::AxisType>();
  if (axisType == Acts::AxisType::Equidistant) {
    const auto nbins = j.at("axis").at("bins").get<std::size_t>();
    if (nbins == 0u) {
      throw std::invalid_argument("Number of bins must be positive");
    }
    if (j.at("autorange").get<bool>()) {
      return Acts::ProtoAxis(axisBoundaryType, nbins);
    }
    const auto min = j.at("axis").at("range").at(0).get<double>();
    const auto max = j.at("axis").at("range").at(1).get<double>();
    if (min >= max) {
      throw std::invalid_argument("Invalid range: min must be less than max");
    }
    return Acts::ProtoAxis(axisBoundaryType, min, max, nbins);
  }
  const auto binEdges =
      j.at("axis").at("boundaries").get<std::vector<double>>();
  if (binEdges.size() < 2u) {
    throw std::invalid_argument("At least two bin edges required");
  }
  if (!std::ranges::is_sorted(binEdges)) {
    throw std::invalid_argument("Bin edges must be sorted in ascending order");
  }
  return Acts::ProtoAxis(axisBoundaryType, binEdges);
}

}  // namespace

void Acts::to_json(nlohmann::json& j, const Acts::BinningData& bd) {
  // Common to all bin utilities
  j["min"] = bd.min;
  j["max"] = bd.max;
  j["option"] = (bd.option == Acts::open ? "open" : "closed");
  j["value"] = bd.binvalue;
  int bins = bd.bins();
  // Write sub bin data if present
  if (bd.subBinningData != nullptr) {
    nlohmann::json subjson;
    to_json(subjson, *bd.subBinningData);
    j["subdata"] = subjson;
    j["subadditive"] = bd.subBinningAdditive;
    // this modifies the bins as bins() returns total number in general
    if (bd.subBinningAdditive) {
      bins -= static_cast<int>(subjson["bins"]) + 1;
    } else {
      bins /= static_cast<int>(subjson["bins"]);
    }
  }
  // Now distinguish between equidistant / arbitrary
  if (bd.type == Acts::equidistant) {
    j["type"] = "equidistant";
  } else if (bd.type == Acts::arbitrary) {
    j["type"] = "arbitrary";
    j["boundaries"] = bd.boundaries();
  }
  j["bins"] = bins;
}

void Acts::from_json(const nlohmann::json& j, BinningData& bd) {
  // Common to all bin utilities
  float min = j["min"];
  float max = j["max"];
  int bins = j["bins"];

  // Support legacy format with BinningValue instead of AxisDirection,
  // this will anyway disappear with the removal of BinUtility
  AxisDirection bValue = AxisDirection::AxisX;
  if (j["value"].get<std::string>().substr(0, 3) == "bin") {
    std::string bValueStr = j["value"];
    if (bValueStr == "binX") {
      bValue = AxisDirection::AxisX;
    } else if (bValueStr == "binY") {
      bValue = AxisDirection::AxisY;
    } else if (bValueStr == "binZ") {
      bValue = AxisDirection::AxisZ;
    } else if (bValueStr == "binR") {
      bValue = AxisDirection::AxisR;
    } else if (bValueStr == "binPhi") {
      bValue = AxisDirection::AxisPhi;
    } else if (bValueStr == "binRPhi") {
      bValue = AxisDirection::AxisRPhi;
    } else if (bValueStr == "binH") {
      bValue = AxisDirection::AxisTheta;
    } else if (bValueStr == "binEta") {
      bValue = AxisDirection::AxisEta;
    } else if (bValueStr == "binMag") {
      bValue = AxisDirection::AxisMag;
    } else {
      throw std::invalid_argument("Unknown binning value name: " + bValueStr);
    }
  } else {
    bValue = j["value"].get<AxisDirection>();
  }

  if (bins == 1 && !(j["type"] == "arbitrary")) {
    bd = BinningData(bValue, min, max);
    return;
  }
  Acts::BinningOption bOption = (j["option"] == "open") ? open : closed;
  Acts::BinningType bType =
      (j["type"] == "equidistant") ? equidistant : arbitrary;

  std::unique_ptr<BinningData> subBinning = nullptr;
  bool subBinningAdditive = false;
  if (j.find("subdata") != j.end()) {
    subBinningAdditive = j["subadditive"];
  }

  if (bType == equidistant) {
    bd = BinningData(bOption, bValue, bins, min, max, std::move(subBinning),
                     subBinningAdditive);
  } else {
    std::vector<float> boundaries = j["boundaries"];
    bd = BinningData(bOption, bValue, boundaries, std::move(subBinning));
  }
}

void Acts::to_json(nlohmann::json& j, const BinUtility& bu) {
  nlohmann::json jbindata;
  for (const auto& bdata : bu.binningData()) {
    jbindata.push_back(nlohmann::json(bdata));
  }
  j["binningdata"] = jbindata;
  if (!bu.transform().isApprox(Transform3::Identity())) {
    nlohmann::json jtrf = Transform3JsonConverter::toJson(bu.transform());
    j["transform"] = jtrf;
  }
}

void Acts::from_json(const nlohmann::json& j, Acts::BinUtility& bu) {
  bu = Acts::BinUtility();
  if (j.find("transform") != j.end() && !j["transform"].empty()) {
    Acts::Transform3 trf = Transform3JsonConverter::fromJson(j["transform"]);
    bu = Acts::BinUtility(trf);
  }
  for (const auto& jdata : j["binningdata"]) {
    Acts::BinningData bd;
    from_json(jdata, bd);
    bu += Acts::BinUtility(bd);
  }
}

void Acts::to_json(nlohmann::json& j, const Acts::ProtoAxis& pa) {
  j["axis"] = AxisJsonConverter::toJson(pa.getAxis());
  j["autorange"] = pa.isAutorange();
}

void Acts::from_json(const nlohmann::json& j, Acts::ProtoAxis& pa) {
  pa = protoAxisFromJsonImpl(j);
}

void Acts::to_json(nlohmann::json& j, const Acts::DirectedProtoAxis& dpa) {
  j["axis"] = AxisJsonConverter::toJson(dpa.getAxis());
  j["autorange"] = dpa.isAutorange();
  j["direction"] = dpa.getAxisDirection();
}

void Acts::from_json(const nlohmann::json& j, Acts::DirectedProtoAxis& dpa) {
  const auto direction = j.at("direction").get<Acts::AxisDirection>();
  const auto axisBoundaryType =
      j.at("axis").at("boundary_type").get<Acts::AxisBoundaryType>();
  const auto axisType = j.at("axis").at("type").get<Acts::AxisType>();
  if (axisType == Acts::AxisType::Equidistant) {
    const auto nbins = j.at("axis").at("bins").get<std::size_t>();
    if (nbins == 0u) {
      throw std::invalid_argument("Number of bins must be positive");
    }
    if (j.at("autorange").get<bool>()) {
      dpa = Acts::DirectedProtoAxis(direction, axisBoundaryType, nbins);
      return;
    }
    const auto min = j.at("axis").at("range").at(0).get<double>();
    const auto max = j.at("axis").at("range").at(1).get<double>();
    if (min >= max) {
      throw std::invalid_argument("Invalid range: min must be less than max");
    }
    dpa = Acts::DirectedProtoAxis(direction, axisBoundaryType, min, max, nbins);
    return;
  }
  const auto binEdges =
      j.at("axis").at("boundaries").get<std::vector<double>>();
  if (binEdges.size() < 2u) {
    throw std::invalid_argument("At least two bin edges required");
  }
  if (!std::ranges::is_sorted(binEdges)) {
    throw std::invalid_argument("Bin edges must be sorted in ascending order");
  }
  dpa = Acts::DirectedProtoAxis(direction, axisBoundaryType, binEdges);
}
