// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Calorimeter/EDM4hepCalorimeterMetadata.hpp"

#include "Acts/Definitions/Units.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <utility>

namespace ActsExamples {
EDM4hepCalorimeterMetadata::EDM4hepCalorimeterMetadata(
    const Config& config, Acts::Logging::Level level)
    : IAlgorithm("EDM4hepCalorimeterMetadata",
                 Acts::getDefaultLogger("EDM4hepCalorimeterMetadata", level)) {
  if (config.inputFrame.empty() || config.outputFrame.empty() ||
      config.inputFrame == config.outputFrame ||
      config.geometryIdentifier.empty() ||
      config.input.inputFrame != config.inputFrame ||
      config.response.inputSimHits != config.input.outputDeposits ||
      config.response.outputHits != config.clustering.inputHits ||
      config.response.outputHits != config.hitOutput.inputHits ||
      config.input.outputDeposits != config.hitOutput.inputDeposits ||
      config.input.outputSources != config.hitOutput.inputSources ||
      config.clustering.inputHits != config.clusterOutput.inputHits ||
      config.clustering.outputClusters != config.clusterOutput.inputClusters ||
      config.hitOutput.outputHits != config.clusterOutput.inputEdmHits) {
    throw std::invalid_argument(
        "Inconsistent calorimeter metadata configuration");
  }
  // Use the same validation as the reconstruction and converters.
  (void)EDM4hepCalorimeterInputConverter{config.input};
  (void)CalorimeterDigitizationAlgorithm{config.response};
  (void)CalorimeterClusteringAlgorithm{config.clustering};
  (void)EDM4hepCalorimeterOutputConverter{config.hitOutput};
  (void)EDM4hepCalorimeterClusterOutputConverter{config.clusterOutput};
  auto& strings = m_parameters.strings;
  auto& numbers = m_parameters.numbers;
  strings["acts.calo.geometry"] = {config.geometryIdentifier};
  strings["acts.calo.units"] = {"energy=GeV", "length=mm", "time=ns"};
  strings["acts.calo.response.algorithm"] = {"CalorimeterResponse"};
  strings["acts.calo.response.conventions"] = {
      "inclusive deposit time window", "inclusive calibrated cell threshold",
      "accepted-energy-weighted time; no TOF correction"};
  strings["acts.calo.clustering.algorithm"] = {"seeded connected components"};
  strings["acts.calo.clustering.conventions"] = {
      "inclusive cell seed/neighbour thresholds",
      "energy sum; no second calibration", "energy-weighted position and time",
      "highest-energy seed; smaller cell ID breaks ties"};
  strings["acts.calo.cluster.unmeasured"] = {"intrinsic direction",
                                             "covariance", "shower shapes"};
  strings["acts.calo.links.weight"] = {
      "accepted deposited energy fraction per calibrated cell"};
  strings["acts.calo.input.simHits"] = {config.input.inputSimHits};
  strings["acts.calo.native.deposits"] = {config.input.outputDeposits};
  strings["acts.calo.native.sources"] = {config.input.outputSources};
  strings["acts.calo.native.hits"] = {config.response.outputHits};
  strings["acts.calo.native.clusters"] = {config.clustering.outputClusters};
  strings["acts.calo.output.hits"] = {config.hitOutput.outputHits};
  strings["acts.calo.output.truthLinks"] = {config.hitOutput.outputLinks};
  strings["acts.calo.output.clusters"] = {config.clusterOutput.outputClusters};
  strings["acts.calo.output.clusterTimes"] = {config.clusterOutput.outputTimes};
  strings["acts.calo.output.clusterSeedCellIds"] = {
      config.clusterOutput.outputSeedCellIds};
  strings["acts.calo.cluster.sidecars"] = {
      "same index as output cluster; time in ns; seed ID uint64"};
  numbers["acts.calo.schemaVersion"] = {1};
  numbers["acts.calo.response.energyScale"] = {
      config.response.response.energyScale};
  numbers["acts.calo.response.energyThreshold"] = {
      config.response.response.energyThreshold / Acts::UnitConstants::GeV};
  numbers["acts.calo.response.timeWindow"] = {
      config.response.response.timeMin / Acts::UnitConstants::ns,
      config.response.response.timeMax / Acts::UnitConstants::ns};
  numbers["acts.calo.clustering.seedEnergyThreshold"] = {
      config.clustering.clustering.seedEnergyThreshold /
      Acts::UnitConstants::GeV};
  numbers["acts.calo.clustering.neighbourEnergyThreshold"] = {
      config.clustering.clustering.neighbourEnergyThreshold /
      Acts::UnitConstants::GeV};
  std::set<std::pair<std::uint64_t, std::uint64_t>> edges;
  for (const auto& edge : config.clustering.neighbours) {
    edges.emplace(std::min(edge.first, edge.second),
                  std::max(edge.first, edge.second));
  }
  auto& topology = strings["acts.calo.clustering.neighbours"];
  for (const auto& [first, second] : edges) {
    topology.push_back(std::to_string(first) + ":" + std::to_string(second));
  }
  std::set<std::string> keys;
  for (const auto& [key, values] : strings) {
    keys.insert(key);
  }
  for (const auto& [key, values] : numbers) {
    keys.insert(key);
  }
  const std::map<std::string, std::string> jetInputs{
      {"acts.calo.jets.inputClusters", config.clustering.outputClusters},
      {"acts.calo.jets.inputEdmClusters", config.clusterOutput.outputClusters},
      {"acts.calo.jets.inputHits", config.response.outputHits},
      {"acts.calo.jets.inputEdmHits", config.hitOutput.outputHits}};
  for (const auto& [key, expected] : jetInputs) {
    const auto found = config.additional.strings.find(key);
    if (config.additional.numbers.contains(key) ||
        (found != config.additional.strings.end() &&
         found->second != std::vector<std::string>{expected})) {
      throw std::invalid_argument("Jet metadata collection names disagree");
    }
  }
  for (const auto& [key, values] : config.additional.strings) {
    if (!key.starts_with("acts.calo.") || !keys.insert(key).second) {
      throw std::invalid_argument(
          "Duplicate or invalid calorimeter metadata key");
    }
    strings.emplace(key, values);
  }
  for (const auto& [key, values] : config.additional.numbers) {
    if (!key.starts_with("acts.calo.") || !keys.insert(key).second ||
        std::ranges::any_of(
            values, [](double value) { return !std::isfinite(value); })) {
      throw std::invalid_argument(
          "Duplicate or invalid calorimeter metadata value");
    }
    numbers.emplace(key, values);
  }
  m_input.initialize(config.inputFrame);
  m_output.initialize(config.outputFrame);
}
ProcessCode EDM4hepCalorimeterMetadata::execute(
    const AlgorithmContext& context) const {
  auto frame = m_input(context);
  auto rejectExisting = [](const auto& keys) {
    if (std::ranges::any_of(keys, [](const auto& key) {
          return key.starts_with("acts.calo.");
        })) {
      throw std::invalid_argument(
          "Input frame already has calorimeter metadata");
    }
  };
  rejectExisting(frame.getParameterKeys<int>());
  rejectExisting(frame.getParameterKeys<float>());
  rejectExisting(frame.getParameterKeys<double>());
  rejectExisting(frame.getParameterKeys<std::string>());
  for (const auto& [key, values] : m_parameters.strings) {
    frame.putParameter(key, values);
  }
  for (const auto& [key, values] : m_parameters.numbers) {
    frame.putParameter(key, values);
  }
  m_output(context, std::move(frame));
  return ProcessCode::SUCCESS;
}
}  // namespace ActsExamples
