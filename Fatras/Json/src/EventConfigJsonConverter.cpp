// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Json/EventConfigJsonConverter.hpp"

#include "Acts/Definitions/PdgParticle.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>

namespace ActsFatras::Synthetic {

namespace {

constexpr const char* kConfigHeader = "acts-fatras-synthetic-config";
constexpr int kFormatVersion = 0;

/// A float as the shortest decimal that reads back as itself; see the same
/// helper in `DetectorDescriptionJsonConverter.cpp` for why a widened float is
/// not what to write.
/// @param value the number to write
/// @return the number
double number(const float value) {
  for (int digits = 1; digits <= 9; ++digits) {
    std::array<char, 32> text{};
    std::snprintf(text.data(), text.size(), "%.*g", digits,
                  static_cast<double>(value));
    const double shortest = std::strtod(text.data(), nullptr);
    if (static_cast<float>(shortest) == value) {
      return shortest;
    }
  }
  return static_cast<double>(value);
}

}  // namespace

void to_json(nlohmann::json& j, const EnergyLossModel& model) {
  j = model == EnergyLossModel::Mean ? "mean" : "mode";
}

void from_json(const nlohmann::json& j, EnergyLossModel& model) {
  const std::string text = j.get<std::string>();
  if (text == "mean") {
    model = EnergyLossModel::Mean;
    return;
  }
  if (text == "mode") {
    model = EnergyLossModel::Mode;
    return;
  }
  throw std::invalid_argument(
      "EventConfig: an energy loss is the 'mean' or the 'mode', not '" + text +
      "'");
}

void to_json(nlohmann::json& j, const GenerationConfig& config) {
  j = nlohmann::json{
      {"pileup", config.pileup},
      {"chargedPerUnitRapidity", number(config.chargedPerUnitRapidity)},
      {"minPt", number(config.minPt)},
      {"ptScale", number(config.ptScale)},
      {"ptExponent", number(config.ptExponent)},
      {"minRapidity", number(config.minRapidity)},
      {"maxRapidity", number(config.maxRapidity)},
      {"rapidityEdge", number(config.rapidityEdge)},
      {"rapidityEdgeWidth", number(config.rapidityEdgeWidth)},
      {"beamspotSigmaZ", number(config.beamspotSigmaZ)},
      {"d0Sigma", number(config.d0Sigma)}};
}

void from_json(const nlohmann::json& j, GenerationConfig& config) {
  config.pileup = j.at("pileup").get<std::size_t>();
  config.chargedPerUnitRapidity = j.at("chargedPerUnitRapidity").get<float>();
  config.minPt = j.at("minPt").get<float>();
  config.ptScale = j.at("ptScale").get<float>();
  config.ptExponent = j.at("ptExponent").get<float>();
  config.minRapidity = j.at("minRapidity").get<float>();
  config.maxRapidity = j.at("maxRapidity").get<float>();
  config.rapidityEdge = j.at("rapidityEdge").get<float>();
  config.rapidityEdgeWidth = j.at("rapidityEdgeWidth").get<float>();
  config.beamspotSigmaZ = j.at("beamspotSigmaZ").get<float>();
  config.d0Sigma = j.at("d0Sigma").get<float>();
}

void to_json(nlohmann::json& j, const PropagationConfig& config) {
  j = nlohmann::json{
      {"maxTurns", number(config.maxTurns)},
      {"tracksPerPrimary", number(config.tracksPerPrimary)},
      {"hitsPerPrimary", number(config.hitsPerPrimary)},
      {"maxTracksPerPrimary", number(config.maxTracksPerPrimary)}};
}

void from_json(const nlohmann::json& j, PropagationConfig& config) {
  config.maxTurns = j.at("maxTurns").get<float>();
  config.tracksPerPrimary = j.at("tracksPerPrimary").get<float>();
  config.hitsPerPrimary = j.at("hitsPerPrimary").get<float>();
  config.maxTracksPerPrimary = j.at("maxTracksPerPrimary").get<float>();
}

void to_json(nlohmann::json& j, const MaterialConfig& config) {
  j = nlohmann::json{
      {"maxDiscPathLength", number(config.maxDiscPathLength)},
      {"maxCylinderPathLength", number(config.maxCylinderPathLength)},
      {"scale", number(config.scale)},
      {"multipleScattering", config.multipleScattering},
      {"energyLoss", config.energyLoss},
      {"energyLossModel", config.energyLossModel},
      {"maxEnergyLossFraction", number(config.maxEnergyLossFraction)}};
}

void from_json(const nlohmann::json& j, MaterialConfig& config) {
  config.maxDiscPathLength = j.at("maxDiscPathLength").get<float>();
  config.maxCylinderPathLength = j.at("maxCylinderPathLength").get<float>();
  config.scale = j.at("scale").get<float>();
  config.multipleScattering = j.at("multipleScattering").get<bool>();
  config.energyLoss = j.at("energyLoss").get<bool>();
  config.energyLossModel = j.at("energyLossModel").get<EnergyLossModel>();
  config.maxEnergyLossFraction = j.at("maxEnergyLossFraction").get<float>();
}

void to_json(nlohmann::json& j, const MeasurementConfig& config) {
  j = nlohmann::json{{"positionSmearing", number(config.positionSmearing)},
                     {"overlapScale", number(config.overlapScale)}};
}

void from_json(const nlohmann::json& j, MeasurementConfig& config) {
  config.positionSmearing = j.at("positionSmearing").get<float>();
  config.overlapScale = j.at("overlapScale").get<float>();
}

void to_json(nlohmann::json& j, const SecondarySamplingConfig& config) {
  j = nlohmann::json{
      {"minPt", number(config.minPt)},
      {"electronScale", number(config.electronScale)},
      {"electronExponent", number(config.electronExponent)},
      {"electronSpread", number(config.electronSpread)},
      {"electronKt", number(config.electronKt)},
      {"momentumScale", number(config.momentumScale)},
      {"momentumExponent", number(config.momentumExponent)},
      {"momentumSpread", number(config.momentumSpread)},
      {"kt", number(config.kt)},
      {"evaporationFraction", number(config.evaporationFraction)},
      {"evaporationScale", number(config.evaporationScale)}};
}

void from_json(const nlohmann::json& j, SecondarySamplingConfig& config) {
  config.minPt = j.at("minPt").get<float>();
  config.electronScale = j.at("electronScale").get<float>();
  config.electronExponent = j.at("electronExponent").get<float>();
  config.electronSpread = j.at("electronSpread").get<float>();
  config.electronKt = j.at("electronKt").get<float>();
  config.momentumScale = j.at("momentumScale").get<float>();
  config.momentumExponent = j.at("momentumExponent").get<float>();
  config.momentumSpread = j.at("momentumSpread").get<float>();
  config.kt = j.at("kt").get<float>();
  config.evaporationFraction = j.at("evaporationFraction").get<float>();
  config.evaporationScale = j.at("evaporationScale").get<float>();
}

void to_json(nlohmann::json& j, const SecondaryConfig& config) {
  j = nlohmann::json{{"electronRate", number(config.electronRate)},
                     {"nuclearRate", number(config.nuclearRate)},
                     {"decayYield", number(config.decayYield)},
                     {"decayLength", number(config.decayLength)},
                     {"stubRate", number(config.stubRate)},
                     {"stubClusters", number(config.stubClusters)},
                     {"stubReach", number(config.stubReach)},
                     {"maxGenerations", config.maxGenerations},
                     {"maxPerCrossing", number(config.maxPerCrossing)},
                     {"sampling", config.sampling}};
}

void from_json(const nlohmann::json& j, SecondaryConfig& config) {
  config.electronRate = j.at("electronRate").get<float>();
  config.nuclearRate = j.at("nuclearRate").get<float>();
  config.decayYield = j.at("decayYield").get<float>();
  config.decayLength = j.at("decayLength").get<float>();
  config.stubRate = j.at("stubRate").get<float>();
  config.stubClusters = j.at("stubClusters").get<float>();
  config.stubReach = j.at("stubReach").get<float>();
  config.maxGenerations = j.at("maxGenerations").get<std::uint32_t>();
  config.maxPerCrossing = j.at("maxPerCrossing").get<float>();
  config.sampling = j.at("sampling").get<SecondarySamplingConfig>();
}

void to_json(nlohmann::json& j, const SimulationConfig& config) {
  j = nlohmann::json{{"propagation", config.propagation},
                     {"material", config.material},
                     {"measurement", config.measurement},
                     {"secondaries", config.secondaries}};
}

void from_json(const nlohmann::json& j, SimulationConfig& config) {
  config.propagation = j.at("propagation").get<PropagationConfig>();
  config.material = j.at("material").get<MaterialConfig>();
  config.measurement = j.at("measurement").get<MeasurementConfig>();
  config.secondaries = j.at("secondaries").get<SecondaryConfig>();
}

void to_json(nlohmann::json& j, const EventConfig& config) {
  j = nlohmann::json{
      {"generation", config.generation},
      {"simulation", config.simulation},
      {"particlePdg", static_cast<std::int32_t>(config.particlePdg)},
      {"bFieldZ", number(config.bFieldZ)},
      {"seed", config.seed}};
}

void from_json(const nlohmann::json& j, EventConfig& config) {
  config.generation = j.at("generation").get<GenerationConfig>();
  config.simulation = j.at("simulation").get<SimulationConfig>();
  config.particlePdg =
      static_cast<Acts::PdgParticle>(j.at("particlePdg").get<std::int32_t>());
  config.bFieldZ = j.at("bFieldZ").get<float>();
  config.seed = j.at("seed").get<std::uint32_t>();
}

}  // namespace ActsFatras::Synthetic

namespace ActsFatras {

Synthetic::EventConfig Synthetic::readEventConfig(
    const std::filesystem::path& path) {
  std::ifstream file(path);
  if (!file.is_open()) {
    throw std::runtime_error("ActsFatras::Synthetic: cannot read " +
                             path.string());
  }
  nlohmann::json j;
  try {
    file >> j;
  } catch (const nlohmann::json::parse_error& error) {
    throw std::runtime_error("ActsFatras::Synthetic: " + path.string() +
                             " is not JSON: " + error.what());
  }
  if (!j.contains(kConfigHeader)) {
    throw std::runtime_error("ActsFatras::Synthetic: " + path.string() +
                             " carries no '" + kConfigHeader +
                             "' header, so it is some other kind of file");
  }
  return j.get<EventConfig>();
}

void Synthetic::writeEventConfig(const std::filesystem::path& path,
                                 const EventConfig& config) {
  nlohmann::json j = config;
  j[kConfigHeader] = nlohmann::json{{"format-version", kFormatVersion}};
  std::ofstream file(path);
  if (!file.is_open()) {
    throw std::runtime_error("ActsFatras::Synthetic: cannot write " +
                             path.string());
  }
  // four spaces, as `CI/format_json.py` normalises to
  file << j.dump(4) << std::endl;
}

}  // namespace ActsFatras
