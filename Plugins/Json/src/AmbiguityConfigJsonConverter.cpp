// This file is part of the Acts project.
//
// Copyright (C) 2023-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/AmbiguityConfigJsonConverter.hpp"

#include "Acts/AmbiguityResolution/ScoreBasedAmbiguityResolution.hpp"
#include "Acts/Plugins/Json/ActsJson.hpp"

namespace Acts {

void from_json(const nlohmann::json& j, ConfigPair& p) {
  std::vector<ScoreBasedAmbiguityResolution::DetectorConfig> detectorConfigs;
  std::map<std::size_t, std::size_t> volumeMap;

  for (auto& [key, value] : j.items()) {
    ScoreBasedAmbiguityResolution::DetectorConfig detectorConfig;

    std::size_t detectorId = std::stoi(key);

    detectorConfig.hitsScoreWeight = value["hitsScoreWeight"];
    detectorConfig.holesScoreWeight = value["holesScoreWeight"];
    detectorConfig.outliersScoreWeight = value["outliersScoreWeight"];
    detectorConfig.otherScoreWeight = value["otherScoreWeight"];

    detectorConfig.minHits = value["minHits"];
    detectorConfig.maxHits = value["maxHits"];
    detectorConfig.maxHoles = value["maxHoles"];
    detectorConfig.maxOutliers = value["maxOutliers"];
    detectorConfig.maxSharedHits = value["maxSharedHits"];

    detectorConfig.sharedHitsFlag = value["sharedHitsFlag"];

    std::vector<double> factorHits = value["factorHits"];
    std::vector<double> factorHoles = value["factorHoles"];

    for (auto factor : factorHits) {
      detectorConfig.factorHits.push_back(factor);
    }

    for (auto factor : factorHoles) {
      detectorConfig.factorHoles.push_back(factor);
    }

    detectorConfigs.push_back(detectorConfig);

    std::vector<std::size_t> volumesIds = value["volumesIds"];
    for (auto volumeId : volumesIds) {
      volumeMap[volumeId] = detectorId;
    }
  }
  p = std::make_pair(volumeMap, detectorConfigs);
}

}  // namespace Acts
