// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/AmbiguityConfigJsonConverter.hpp"

#include "Acts/AmbiguityResolution/ScoreBasedAmbiguityResolution.hpp"
#include "Acts/Plugins/Json/ActsJson.hpp"

void initializeEtaVector(std::vector<std::size_t>& target,
                      const std::vector<std::size_t>& source) {
  target = {};
  for (auto value : source) {
    target.push_back(value);
  }
}
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

    detectorConfig.sharedHitsFlag = value["sharedHitsFlag"];

    const std::vector<double>& goodHits = value["goodHits"];
    const std::vector<double>& goodHoles = value["goodHoles"];

    const std::vector<double>& fakeHits = value["fakeHits"];
    const std::vector<double>& fakeHoles = value["fakeHoles"];

    const std::vector<double>& etaBins = value["etaBins"];

    // Validate eta bins
    if (!etaBins.empty()) {
      // Verify monotonically increasing eta bins
      if (!std::is_sorted(etaBins.begin(), etaBins.end())) {
        throw std::invalid_argument(
            "Eta bins must be monotonically increasing");
      }

      detectorConfig.etaBins = {};
      for (auto etaBin : etaBins) {
        detectorConfig.etaBins.push_back(etaBin);
      }
    }

    const std::vector<std::size_t>& minHitsPerEta = value["minHitsPerEta"];
    initializeEtaVector(detectorConfig.minHitsPerEta, minHitsPerEta);

    const std::vector<std::size_t>& maxHolesPerEta = value["maxHolesPerEta"];
    initializeEtaVector(detectorConfig.maxHolesPerEta, maxHolesPerEta);

    const std::vector<std::size_t>& maxOutliersPerEta =
        value["maxOutliersPerEta"];
    initializeEtaVector(detectorConfig.maxOutliersPerEta, maxOutliersPerEta);

    const std::vector<std::size_t>& maxSharedHitsPerEta =
        value["maxSharedHitsPerEta"];
    initializeEtaVector(detectorConfig.maxSharedHitsPerEta, maxSharedHitsPerEta);

    if (goodHits.size() != fakeHits.size()) {
      throw std::invalid_argument("goodHits and FakeHits size mismatch");
    }

    detectorConfig.factorHits = {};
    detectorConfig.maxHits = goodHits.size() - 1;
    for (std::size_t i = 0; i < goodHits.size(); i++) {
      detectorConfig.factorHits.push_back(goodHits[i] / fakeHits[i]);
    }

    if (goodHoles.size() != fakeHoles.size()) {
      throw std::invalid_argument("goodHoles and FakeHoles size mismatch");
    }

    detectorConfig.factorHoles = {};
    detectorConfig.maxHoles = goodHoles.size() - 1;
    for (std::size_t i = 0; i < goodHoles.size(); i++) {
      detectorConfig.factorHoles.push_back(goodHoles[i] / fakeHoles[i]);
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
