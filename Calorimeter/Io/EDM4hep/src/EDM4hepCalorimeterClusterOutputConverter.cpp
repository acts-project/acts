// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Calorimeter/EDM4hepCalorimeterClusterOutputConverter.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <utility>

#include "EDM4hepCalorimeterValidation.hpp"

namespace ActsExamples {
EDM4hepCalorimeterClusterOutputConverter::
    EDM4hepCalorimeterClusterOutputConverter(const Config& config,
                                             Acts::Logging::Level level)
    : PodioOutputConverter(
          "EDM4hepCalorimeterClusterOutput",
          Acts::getDefaultLogger("EDM4hepCalorimeterClusterOutput", level)),
      m_config(config) {
  const std::set<std::string> names{
      config.inputClusters,  config.inputHits,   config.inputEdmHits,
      config.outputClusters, config.outputTimes, config.outputSeedCellIds};
  if (names.size() != 6 || names.contains("")) {
    throw std::invalid_argument(
        "Cluster converter names must be distinct and nonempty");
  }
  m_clusters.initialize(config.inputClusters);
  m_hits.initialize(config.inputHits);
  m_edmHits.initialize(config.inputEdmHits);
  m_output.initialize(config.outputClusters);
  m_times.initialize(config.outputTimes);
  m_seeds.initialize(config.outputSeedCellIds);
}
ProcessCode EDM4hepCalorimeterClusterOutputConverter::execute(
    const AlgorithmContext& context) const {
  const auto& clusters = m_clusters(context);
  const auto& hits = m_hits(context);
  const auto& edmHits = m_edmHits(context);
  detail::checkCalorimeterHitMapping(hits, edmHits);
  edm4hep::ClusterCollection output;
  podio::UserDataCollection<double> times;
  podio::UserDataCollection<std::uint64_t> seeds;
  std::set<std::size_t> usedHits;
  std::set<std::uint64_t> cells;
  for (const auto& hit : hits) {
    if (!cells.insert(hit.cellId).second) {
      throw std::invalid_argument("Duplicate native calorimeter cell ID");
    }
  }
  for (const auto& cluster : clusters) {
    if (!std::isfinite(cluster.energy) || cluster.energy <= 0 ||
        !cluster.position.allFinite() || !std::isfinite(cluster.time) ||
        cluster.hitIndices.empty()) {
      throw std::invalid_argument("Invalid native calorimeter cluster");
    }
    double energy = 0;
    double time = 0;
    Acts::Vector3 position = Acts::Vector3::Zero();
    const ActsCalorimeter::CalorimeterHit* seed = nullptr;
    double positionScale = 0;
    double timeScale = 0;
    for (const auto index : cluster.hitIndices) {
      if (index >= hits.size() || !usedHits.insert(index).second) {
        throw std::invalid_argument("Invalid or repeated cluster hit index");
      }
      const auto& hit = hits[index];
      const double total = energy + hit.energy;
      if (!std::isfinite(total)) {
        throw std::overflow_error("Cluster constituent energy sum overflow");
      }
      const double fraction = hit.energy / total;
      position = (1 - fraction) * position + fraction * hit.position;
      time = (1 - fraction) * time + fraction * hit.time;
      energy = total;
      positionScale =
          std::max(positionScale, hit.position.cwiseAbs().maxCoeff());
      timeScale = std::max(timeScale, std::abs(hit.time));
      if (seed == nullptr || hit.energy > seed->energy ||
          (hit.energy == seed->energy && hit.cellId < seed->cellId)) {
        seed = &hit;
      }
    }
    if (!detail::calorimeterClose(cluster.energy, energy, energy) ||
        !detail::calorimeterClose(cluster.time, time, timeScale) ||
        !(cluster.position - position).allFinite() ||
        (cluster.position - position).cwiseAbs().maxCoeff() >
            64 * std::numeric_limits<double>::epsilon() * positionScale ||
        cluster.seedCellId != seed->cellId) {
      throw std::invalid_argument(
          "Cluster values disagree with constituent hits");
    }
    auto result = output.create();
    result.setEnergy(detail::calorimeterEnergy(cluster.energy));
    result.setPosition({detail::calorimeterFloat(cluster.position.x() /
                                                 Acts::UnitConstants::mm),
                        detail::calorimeterFloat(cluster.position.y() /
                                                 Acts::UnitConstants::mm),
                        detail::calorimeterFloat(cluster.position.z() /
                                                 Acts::UnitConstants::mm)});
    // Intrinsic shower directions, covariance and shape parameters are not
    // measured by this prototype and retain the EDM defaults.
    for (const auto index : cluster.hitIndices) {
      result.addToHits(edmHits[index]);
    }
    times.push_back(cluster.time / Acts::UnitConstants::ns);
    seeds.push_back(cluster.seedCellId);
  }
  m_output(context, std::move(output));
  m_times(context, std::move(times));
  m_seeds(context, std::move(seeds));
  return ProcessCode::SUCCESS;
}
std::vector<std::string> EDM4hepCalorimeterClusterOutputConverter::collections()
    const {
  return {m_config.outputClusters, m_config.outputTimes,
          m_config.outputSeedCellIds};
}
}  // namespace ActsExamples
