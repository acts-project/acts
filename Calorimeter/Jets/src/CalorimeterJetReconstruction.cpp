// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsCalorimeter/Jets/CalorimeterJetReconstruction.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include <fastjet/ClusterSequence.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

namespace ActsCalorimeter {

CalorimeterJetReconstruction::CalorimeterJetReconstruction(const Config& config)
    : m_config(config) {
  if (!std::isfinite(config.radius) || config.radius <= 0 ||
      config.radius > fastjet::JetDefinition::max_allowable_R ||
      !std::isfinite(config.jetPtMin) || config.jetPtMin < 0 ||
      !config.origin.allFinite()) {
    throw std::invalid_argument("Invalid calorimeter jet configuration");
  }
}

CalorimeterJetContainer CalorimeterJetReconstruction::operator()(
    std::span<const CalorimeterCluster> clusters) const {
  if (clusters.size() >
      static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::overflow_error(
        "Cluster indices exceed FastJet's user-index range");
  }

  std::vector<fastjet::PseudoJet> inputs;
  inputs.reserve(clusters.size());
  double totalEnergy = 0;
  // FastJet computes squared momenta. Keep even recombined momenta comfortably
  // within that range, rather than accepting finite energies that square to
  // inf.
  const double energyLimit = std::sqrt(std::numeric_limits<double>::max()) / 2;
  for (std::size_t index = 0; index < clusters.size(); ++index) {
    const auto& cluster = clusters[index];
    if (!std::isfinite(cluster.energy) || cluster.energy < 0 ||
        !cluster.position.allFinite()) {
      throw std::invalid_argument("Invalid calorimeter cluster kinematics");
    }
    if (cluster.energy == 0) {
      continue;
    }
    const Acts::Vector3 displacement = cluster.position - m_config.origin;
    const double distance = displacement.stableNorm();
    if (!displacement.allFinite() || !std::isfinite(distance) ||
        distance == 0) {
      throw std::invalid_argument("Undefined calorimeter cluster direction");
    }
    totalEnergy += cluster.energy;
    if (!std::isfinite(totalEnergy) || totalEnergy > energyLimit) {
      throw std::overflow_error(
          "Calorimeter energy exceeds FastJet's numeric range");
    }
    const Acts::Vector3 momentum = cluster.energy * (displacement / distance);
    auto& input = inputs.emplace_back(momentum.x(), momentum.y(), momentum.z(),
                                      cluster.energy);
    input.set_user_index(static_cast<int>(index));
  }
  if (inputs.empty()) {
    return {};
  }

  const fastjet::JetDefinition definition(fastjet::antikt_algorithm,
                                          m_config.radius, fastjet::E_scheme);
  const fastjet::ClusterSequence sequence(inputs, definition);
  CalorimeterJetContainer output;
  // Apply the cut with hypot below: inclusive_jets(ptmin) squares its argument.
  for (const auto& jet : sequence.inclusive_jets()) {
    CalorimeterJet result;
    result.fourMomentum = Acts::Vector4{jet.px(), jet.py(), jet.pz(), jet.E()};
    if (!result.fourMomentum.allFinite()) {
      throw std::overflow_error("Nonfinite calorimeter jet four-momentum");
    }
    if (std::hypot(jet.px(), jet.py()) < m_config.jetPtMin) {
      continue;
    }
    for (const auto& constituent : jet.constituents()) {
      const int index = constituent.user_index();
      if (index < 0 || static_cast<std::size_t>(index) >= clusters.size()) {
        throw std::logic_error("Invalid FastJet calorimeter constituent index");
      }
      result.clusterIndices.push_back(static_cast<std::size_t>(index));
    }
    std::ranges::sort(result.clusterIndices);
    output.push_back(std::move(result));
  }
  std::ranges::sort(output, [](const auto& first, const auto& second) {
    const auto firstPt =
        std::hypot(first.fourMomentum[0], first.fourMomentum[1]);
    const auto secondPt =
        std::hypot(second.fourMomentum[0], second.fourMomentum[1]);
    if (firstPt != secondPt) {
      return firstPt > secondPt;
    }
    return first.clusterIndices < second.clusterIndices;
  });
  return output;
}

}  // namespace ActsCalorimeter
