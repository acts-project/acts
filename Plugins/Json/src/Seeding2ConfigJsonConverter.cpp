// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/Seeding2ConfigJsonConverter.hpp"

#include "Acts/Plugins/Json/DefinitionsJsonConverter.hpp"
#include "Acts/Utilities/GridBinFinder.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <algorithm>

namespace {

template <std::size_t DIM>
void to_json(nlohmann::json& j, const Acts::GridBinFinder<DIM>& f) {
  nlohmann::json::array_t values;

  for (const auto& value : f.values()) {
    std::visit(
        Acts::overloaded{
            [&](int v) { values.push_back(v); },
            [&](const std::pair<int, int>& p) { values.push_back(p); },
            [&](const std::vector<std::pair<int, int>>& vec) {
              values.push_back(vec);
            },
            [](const auto&) {
              throw std::runtime_error("Unsupported type in GridBinFinder");
            }},
        value);
  }

  j = values;
}

template <std::size_t DIM>
Acts::GridBinFinder<DIM> from_json(const nlohmann::json& j) {
  using stored_values_t = typename Acts::GridBinFinder<DIM>::stored_values_t;

  std::array<stored_values_t, DIM> values;
  for (std::size_t i = 0; i < DIM; ++i) {
    if (j[i].is_number_integer()) {
      values[i] = j[i].get<int>();
    } else if (j[i].is_array() && j[i].size() == 2 &&
               j[i][0].is_number_integer() && j[i][1].is_number_integer()) {
      values[i] = j[i].get<std::pair<int, int>>();
    } else if (j[i].is_array() &&
               std::ranges::all_of(j[i], [](const nlohmann::json& elem) {
                 return elem.is_array() && elem.size() == 2 &&
                        elem[0].is_number_integer() &&
                        elem[1].is_number_integer();
               })) {
      values[i] = j[i].get<std::vector<std::pair<int, int>>>();
    } else {
      throw std::runtime_error("Invalid type for GridBinFinder value");
    }
  }

  return Acts::GridBinFinder<DIM>(std::move(values));
}

}  // namespace

void Acts::to_json(nlohmann::json& j,
                   const SeedConfirmationRangeConfig& config) {
  j["zMinSeedConf"] = config.zMinSeedConf;
  j["zMaxSeedConf"] = config.zMaxSeedConf;
  j["rMaxSeedConf"] = config.rMaxSeedConf;
  j["nTopForLargeR"] = config.nTopForLargeR;
  j["nTopForSmallR"] = config.nTopForSmallR;
  j["seedConfMinBottomRadius"] = config.seedConfMinBottomRadius;
  j["seedConfMaxZOrigin"] = config.seedConfMaxZOrigin;
  j["minImpactSeedConf"] = config.minImpactSeedConf;
}

void Acts::from_json(const nlohmann::json& j,
                     SeedConfirmationRangeConfig& config) {
  j["zMinSeedConf"].get_to(config.zMinSeedConf);
  j["zMaxSeedConf"].get_to(config.zMaxSeedConf);
  j["rMaxSeedConf"].get_to(config.rMaxSeedConf);
  j["nTopForLargeR"].get_to(config.nTopForLargeR);
  j["nTopForSmallR"].get_to(config.nTopForSmallR);
  j["seedConfMinBottomRadius"].get_to(config.seedConfMinBottomRadius);
  j["seedConfMaxZOrigin"].get_to(config.seedConfMaxZOrigin);
  j["minImpactSeedConf"].get_to(config.minImpactSeedConf);
}

void Acts::Experimental::to_json(nlohmann::json& j,
                                 const DoubletSeedFinder::Config& config) {
  j["candidateDirection"] = config.candidateDirection;
  j["deltaRMin"] = config.deltaRMin;
  j["deltaRMax"] = config.deltaRMax;
  j["deltaZMin"] = config.deltaZMin;
  j["deltaZMax"] = config.deltaZMax;
  j["impactMax"] = config.impactMax;
  j["interactionPointCut"] = config.interactionPointCut;
  j["collisionRegionMin"] = config.collisionRegionMin;
  j["collisionRegionMax"] = config.collisionRegionMax;
  j["cotThetaMax"] = config.cotThetaMax;
  j["minPt"] = config.minPt;
  j["helixCutTolerance"] = config.helixCutTolerance;
  // experiment cuts cannot be serialized directly, so we skip it
}

void Acts::Experimental::to_json(
    nlohmann::json& j, const DoubletSeedFinder::DerivedConfig& config) {
  to_json(j, static_cast<const DoubletSeedFinder::Config&>(config));
  j["bFieldInZ"] = config.bFieldInZ;
  j["minHelixDiameter2"] = config.minHelixDiameter2;
}

void Acts::Experimental::to_json(
    nlohmann::json& j, const BroadTripletSeedFinder::Options& options) {
  j["bFieldInZ"] = options.bFieldInZ;
  j["useStripMeasurementInfo"] = options.useStripMeasurementInfo;
}

void Acts::Experimental::to_json(
    nlohmann::json& j, const BroadTripletSeedFinder::TripletCuts& cuts) {
  j["minPt"] = cuts.minPt;
  j["sigmaScattering"] = cuts.sigmaScattering;
  j["radLengthPerSeed"] = cuts.radLengthPerSeed;
  j["maxPtScattering"] = cuts.maxPtScattering;
  j["impactMax"] = cuts.impactMax;
  j["helixCutTolerance"] = cuts.helixCutTolerance;
  j["toleranceParam"] = cuts.toleranceParam;
}

void Acts::Experimental::to_json(
    nlohmann::json& j, const BroadTripletSeedFinder::DerivedTripletCuts& cuts) {
  to_json(j, static_cast<const BroadTripletSeedFinder::TripletCuts&>(cuts));
  j["bFieldInZ"] = cuts.bFieldInZ;
  j["highland"] = cuts.highland;
  j["pTPerHelixRadius"] = cuts.pTPerHelixRadius;
  j["minHelixDiameter2"] = cuts.minHelixDiameter2;
  j["sigmapT2perRadius"] = cuts.sigmapT2perRadius;
  j["multipleScattering2"] = cuts.multipleScattering2;
}

void Acts::Experimental::to_json(nlohmann::json& j,
                                 const BroadTripletSeedFilter::Config& config) {
  j["deltaInvHelixDiameter"] = config.deltaInvHelixDiameter;
  j["deltaRMin"] = config.deltaRMin;
  j["compatSeedWeight"] = config.compatSeedWeight;
  j["impactWeightFactor"] = config.impactWeightFactor;
  j["zOriginWeightFactor"] = config.zOriginWeightFactor;
  j["maxSeedsPerSpM"] = config.maxSeedsPerSpM;
  j["compatSeedLimit"] = config.compatSeedLimit;
  j["seedWeightIncrement"] = config.seedWeightIncrement;
  j["numSeedIncrement"] = config.numSeedIncrement;
  j["seedConfirmation"] = config.seedConfirmation;
  j["centralSeedConfirmationRange"] = config.centralSeedConfirmationRange;
  j["forwardSeedConfirmationRange"] = config.forwardSeedConfirmationRange;
  j["maxSeedsPerSpMConf"] = config.maxSeedsPerSpMConf;
  j["maxQualitySeedsPerSpMConf"] = config.maxQualitySeedsPerSpMConf;
  j["useDeltaRinsteadOfTopRadius"] = config.useDeltaRinsteadOfTopRadius;
  // experiment cuts cannot be serialized directly, so we skip it
}

void Acts::Experimental::to_json(
    nlohmann::json& j, const CylindricalSpacePointGrid2::Config& config) {
  j["minPt"] = config.minPt;
  j["rMin"] = config.rMin;
  j["rMax"] = config.rMax;
  j["zMin"] = config.zMin;
  j["zMax"] = config.zMax;
  j["deltaRMax"] = config.deltaRMax;
  j["cotThetaMax"] = config.cotThetaMax;
  j["impactMax"] = config.impactMax;
  j["phiMin"] = config.phiMin;
  j["phiMax"] = config.phiMax;
  j["phiBinDeflectionCoverage"] = config.phiBinDeflectionCoverage;
  j["maxPhiBins"] = config.maxPhiBins;
  j["zBinEdges"] = config.zBinEdges;
  j["rBinEdges"] = config.rBinEdges;
  j["bFieldInZ"] = config.bFieldInZ;
  if (config.bottomBinFinder.has_value()) {
    ::to_json(j["bottomBinFinder"], config.bottomBinFinder.value());
  }
  if (config.topBinFinder.has_value()) {
    ::to_json(j["topBinFinder"], config.topBinFinder.value());
  }
  j["navigation"] = config.navigation;
}

void Acts::Experimental::from_json(const nlohmann::json& j,
                                   DoubletSeedFinder::Config& config) {
  j["candidateDirection"].get_to(config.candidateDirection);
  j["deltaRMin"].get_to(config.deltaRMin);
  j["deltaRMax"].get_to(config.deltaRMax);
  j["deltaZMin"].get_to(config.deltaZMin);
  j["deltaZMax"].get_to(config.deltaZMax);
  j["impactMax"].get_to(config.impactMax);
  j["interactionPointCut"].get_to(config.interactionPointCut);
  j["collisionRegionMin"].get_to(config.collisionRegionMin);
  j["collisionRegionMax"].get_to(config.collisionRegionMax);
  j["cotThetaMax"].get_to(config.cotThetaMax);
  j["minPt"].get_to(config.minPt);
  j["helixCutTolerance"].get_to(config.helixCutTolerance);
  // experiment cuts cannot be serialized directly, so we skip it
}

void Acts::Experimental::from_json(const nlohmann::json& j,
                                   DoubletSeedFinder::DerivedConfig& config) {
  from_json(j, static_cast<DoubletSeedFinder::Config&>(config));
  j["bFieldInZ"].get_to(config.bFieldInZ);
  j["minHelixDiameter2"].get_to(config.minHelixDiameter2);
}

void Acts::Experimental::from_json(const nlohmann::json& j,
                                   BroadTripletSeedFinder::Options& options) {
  j["bFieldInZ"].get_to(options.bFieldInZ);
  j["useStripMeasurementInfo"].get_to(options.useStripMeasurementInfo);
}

void Acts::Experimental::from_json(const nlohmann::json& j,
                                   BroadTripletSeedFinder::TripletCuts& cuts) {
  j["minPt"].get_to(cuts.minPt);
  j["sigmaScattering"].get_to(cuts.sigmaScattering);
  j["radLengthPerSeed"].get_to(cuts.radLengthPerSeed);
  j["maxPtScattering"].get_to(cuts.maxPtScattering);
  j["impactMax"].get_to(cuts.impactMax);
  j["helixCutTolerance"].get_to(cuts.helixCutTolerance);
  j["toleranceParam"].get_to(cuts.toleranceParam);
}

void Acts::Experimental::from_json(
    const nlohmann::json& j, BroadTripletSeedFinder::DerivedTripletCuts& cuts) {
  from_json(j, static_cast<BroadTripletSeedFinder::TripletCuts&>(cuts));
  j["bFieldInZ"].get_to(cuts.bFieldInZ);
  j["highland"].get_to(cuts.highland);
  j["pTPerHelixRadius"].get_to(cuts.pTPerHelixRadius);
  j["minHelixDiameter2"].get_to(cuts.minHelixDiameter2);
  j["sigmapT2perRadius"].get_to(cuts.sigmapT2perRadius);
  j["multipleScattering2"].get_to(cuts.multipleScattering2);
}

void Acts::Experimental::from_json(const nlohmann::json& j,
                                   BroadTripletSeedFilter::Config& config) {
  j["deltaInvHelixDiameter"].get_to(config.deltaInvHelixDiameter);
  j["deltaRMin"].get_to(config.deltaRMin);
  j["compatSeedWeight"].get_to(config.compatSeedWeight);
  j["impactWeightFactor"].get_to(config.impactWeightFactor);
  j["zOriginWeightFactor"].get_to(config.zOriginWeightFactor);
  j["maxSeedsPerSpM"].get_to(config.maxSeedsPerSpM);
  j["compatSeedLimit"].get_to(config.compatSeedLimit);
  j["seedWeightIncrement"].get_to(config.seedWeightIncrement);
  j["numSeedIncrement"].get_to(config.numSeedIncrement);
  j["seedConfirmation"].get_to(config.seedConfirmation);
  j["centralSeedConfirmationRange"].get_to(config.centralSeedConfirmationRange);
  j["forwardSeedConfirmationRange"].get_to(config.forwardSeedConfirmationRange);
  j["maxSeedsPerSpMConf"].get_to(config.maxSeedsPerSpMConf);
  j["maxQualitySeedsPerSpMConf"].get_to(config.maxQualitySeedsPerSpMConf);
  j["useDeltaRinsteadOfTopRadius"].get_to(config.useDeltaRinsteadOfTopRadius);
  // experiment cuts cannot be serialized directly, so we skip it
}

void Acts::Experimental::from_json(const nlohmann::json& j,
                                   CylindricalSpacePointGrid2::Config& config) {
  j["minPt"].get_to(config.minPt);
  j["rMin"].get_to(config.rMin);
  j["rMax"].get_to(config.rMax);
  j["zMin"].get_to(config.zMin);
  j["zMax"].get_to(config.zMax);
  j["deltaRMax"].get_to(config.deltaRMax);
  j["cotThetaMax"].get_to(config.cotThetaMax);
  j["impactMax"].get_to(config.impactMax);
  j["phiMin"].get_to(config.phiMin);
  j["phiMax"].get_to(config.phiMax);
  j["phiBinDeflectionCoverage"].get_to(config.phiBinDeflectionCoverage);
  j["maxPhiBins"].get_to(config.maxPhiBins);
  j["zBinEdges"].get_to(config.zBinEdges);
  j["rBinEdges"].get_to(config.rBinEdges);
  j["bFieldInZ"].get_to(config.bFieldInZ);
  if (j.contains("bottomBinFinder")) {
    config.bottomBinFinder = ::from_json<3ul>(j["bottomBinFinder"]);
  }
  if (j.contains("topBinFinder")) {
    config.topBinFinder = ::from_json<3ul>(j["topBinFinder"]);
  }
  j["navigation"].get_to(config.navigation);
}
