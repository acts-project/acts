// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/ProtoDetectorJsonConverter.hpp"

#include "Acts/Plugins/Json/ExtentJsonConverter.hpp"
#include "Acts/Plugins/Json/UtilitiesJsonConverter.hpp"
#include "Acts/Surfaces/Surface.hpp"

void Acts::to_json(nlohmann::json& j, const Acts::ProtoVolume& pv) {
  j["name"] = pv.name;
  j["extent"] = pv.extent;
  j["layerContainer"] = pv.layerContainer;
  j["layerType"] = static_cast<int>(pv.layerType);

  // Helper m ethod to write binnings
  auto writeBinning = [&](const std::vector<BinningData>& binning,
                          const std::string& key) -> void {
    nlohmann::json jbinning;
    for (const auto& bd : binning) {
      jbinning.push_back(bd);
    }
    j[key] = jbinning;
  };
  writeBinning(pv.layerSurfaceBinning, "layerSurfaceBinning");
  nlohmann::json constituents;
  for (const auto& pvc : pv.constituentVolumes) {
    constituents.push_back(pvc);
  }
  j["constituents"] = constituents;
  writeBinning(pv.constituentBinning, "constituentBinning");
}

void Acts::from_json(const nlohmann::json& j, Acts::ProtoVolume& pv) {
  pv.name = j["name"];
  pv.extent = j["extent"];
  pv.layerContainer = j["layerContainer"];
  pv.layerType = static_cast<Surface::SurfaceType>(j["layerType"]);

  // Helper method to read binnings
  auto readBinning = [&](std::vector<BinningData>& binning,
                         const std::string& key) -> void {
    for (const auto& jbinning : j[key]) {
      binning.push_back(jbinning);
    }
  };
  readBinning(pv.layerSurfaceBinning, "layerSurfaceBinning");
  for (const auto& jc : j["constituents"]) {
    pv.constituentVolumes.push_back(jc);
  }
  readBinning(pv.constituentBinning, "constituentBinning");
}

void Acts::to_json(nlohmann::json& j, const Acts::ProtoDetector& pd) {
  j["name"] = pd.name;
  j["world"] = pd.worldVolume;
}

void Acts::from_json(const nlohmann::json& j, Acts::ProtoDetector& pd) {
  pd.name = j["name"];
  pd.worldVolume = j["world"];
}
