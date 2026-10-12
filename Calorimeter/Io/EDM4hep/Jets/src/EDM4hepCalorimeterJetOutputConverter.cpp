// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Calorimeter/EDM4hepCalorimeterJetOutputConverter.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
#include <utility>

#include "EDM4hepCalorimeterValidation.hpp"

namespace ActsExamples {
EDM4hepCalorimeterJetOutputConverter::EDM4hepCalorimeterJetOutputConverter(
    const Config& config, Acts::Logging::Level level)
    : PodioOutputConverter(
          "EDM4hepCalorimeterJetOutput",
          Acts::getDefaultLogger("EDM4hepCalorimeterJetOutput", level)),
      m_config(config) {
  const std::set<std::string> names{
      config.inputJets,    config.inputClusters,    config.inputHits,
      config.inputEdmHits, config.inputEdmClusters, config.outputJets};
  if (names.size() != 6 || names.contains("")) {
    throw std::invalid_argument(
        "Jet converter names must be distinct and nonempty");
  }
  (void)ActsCalorimeter::CalorimeterJetReconstruction{config.reconstruction};
  m_jets.initialize(config.inputJets);
  m_clusters.initialize(config.inputClusters);
  m_hits.initialize(config.inputHits);
  m_edmHits.initialize(config.inputEdmHits);
  m_edmClusters.initialize(config.inputEdmClusters);
  m_output.initialize(config.outputJets);
}
ProcessCode EDM4hepCalorimeterJetOutputConverter::execute(
    const AlgorithmContext& context) const {
  const auto& jets = m_jets(context);
  const auto& clusters = m_clusters(context);
  const auto& hits = m_hits(context);
  const auto& edmHits = m_edmHits(context);
  const auto& edmClusters = m_edmClusters(context);
  detail::checkCalorimeterHitMapping(hits, edmHits);
  if (clusters.size() != edmClusters.size()) {
    throw std::invalid_argument(
        "Native/persistent cluster collection sizes disagree");
  }
  for (std::size_t index = 0; index < clusters.size(); ++index) {
    const auto& cluster = clusters[index];
    const auto output = edmClusters[index];
    if (!std::isfinite(cluster.energy) || cluster.energy <= 0 ||
        !cluster.position.allFinite() || !output.isAvailable() ||
        output.getEnergy() != detail::calorimeterEnergy(cluster.energy) ||
        output.getPosition().x !=
            detail::calorimeterFloat(cluster.position.x() /
                                     Acts::UnitConstants::mm) ||
        output.getPosition().y !=
            detail::calorimeterFloat(cluster.position.y() /
                                     Acts::UnitConstants::mm) ||
        output.getPosition().z !=
            detail::calorimeterFloat(cluster.position.z() /
                                     Acts::UnitConstants::mm) ||
        output.hits_size() != cluster.hitIndices.size()) {
      throw std::invalid_argument(
          "Native/persistent cluster mapping is inconsistent");
    }
    for (std::size_t member = 0; member < cluster.hitIndices.size(); ++member) {
      const auto hitIndex = cluster.hitIndices[member];
      if (hitIndex >= hits.size() ||
          output.getHits(member) != edmHits[hitIndex]) {
        throw std::invalid_argument(
            "Persistent cluster points to the wrong hit");
      }
    }
  }
  edm4hep::ReconstructedParticleCollection output;
  std::set<std::size_t> usedClusters;
  for (const auto& jet : jets) {
    if (!jet.fourMomentum.allFinite() || jet.fourMomentum[3] <= 0 ||
        jet.clusterIndices.empty()) {
      throw std::invalid_argument("Invalid native calorimeter jet");
    }
    Acts::Vector4 expected = Acts::Vector4::Zero();
    for (const auto index : jet.clusterIndices) {
      if (index >= clusters.size() || !usedClusters.insert(index).second) {
        throw std::invalid_argument("Invalid or repeated jet cluster index");
      }
      const auto& cluster = clusters[index];
      const Acts::Vector3 displacement =
          cluster.position - m_config.reconstruction.origin;
      const double distance = displacement.stableNorm();
      if (!displacement.allFinite() || !std::isfinite(distance) ||
          distance == 0) {
        throw std::invalid_argument("Undefined jet constituent direction");
      }
      expected.head<3>() += cluster.energy * (displacement / distance);
      expected[3] += cluster.energy;
    }
    if (!expected.allFinite()) {
      throw std::overflow_error("Jet constituent sum overflow");
    }
    // FastJet and this check can sum constituents in different orders.
    const double tolerance = 64 * std::numeric_limits<double>::epsilon() *
                             static_cast<double>(jet.clusterIndices.size());
    if (!std::isfinite((jet.fourMomentum - expected).cwiseAbs().maxCoeff()) ||
        (jet.fourMomentum - expected).cwiseAbs().maxCoeff() >
            tolerance * expected[3] ||
        std::hypot(jet.fourMomentum[0], jet.fourMomentum[1]) <
            m_config.reconstruction.jetPtMin) {
      throw std::invalid_argument(
          "Jet four-momentum disagrees with E-scheme constituents");
    }
    const double energy = jet.fourMomentum[3];
    const double massSquaredFraction =
        1 - (jet.fourMomentum.head<3>() / energy).squaredNorm();
    if (massSquaredFraction < -tolerance) {
      throw std::invalid_argument("Spacelike calorimeter jet");
    }
    auto result = output.create();
    result.setEnergy(detail::calorimeterEnergy(energy));
    result.setMomentum({detail::calorimeterFloat(jet.fourMomentum[0] /
                                                 Acts::UnitConstants::GeV),
                        detail::calorimeterFloat(jet.fourMomentum[1] /
                                                 Acts::UnitConstants::GeV),
                        detail::calorimeterFloat(jet.fourMomentum[2] /
                                                 Acts::UnitConstants::GeV)});
    result.setMass(detail::calorimeterFloat(
        energy * std::sqrt(std::max(0., massSquaredFraction)) /
        Acts::UnitConstants::GeV));
    for (const auto index : jet.clusterIndices) {
      result.addToClusters(edmClusters[index]);
    }
  }
  m_output(context, std::move(output));
  return ProcessCode::SUCCESS;
}
std::vector<std::string> EDM4hepCalorimeterJetOutputConverter::collections()
    const {
  return {m_config.outputJets};
}
EDM4hepCalorimeterMetadata::Parameters
EDM4hepCalorimeterJetOutputConverter::metadata() const {
  EDM4hepCalorimeterMetadata::Parameters result;
  result.strings["acts.calo.jets.algorithm"] = {"FastJet inclusive anti-kt"};
  result.strings["acts.calo.jets.conventions"] = {
      "massless cluster inputs from fixed origin", "rapidity-phi radius",
      "E-scheme four-vector sum", "inclusive minimum pT",
      "descending pT; constituent indices break ties"};
  result.strings["acts.calo.jets.unmeasured"] = {"PDG", "charge", "covariance"};
  result.strings["acts.calo.jets.inputClusters"] = {m_config.inputClusters};
  result.strings["acts.calo.jets.inputEdmClusters"] = {
      m_config.inputEdmClusters};
  result.strings["acts.calo.jets.inputHits"] = {m_config.inputHits};
  result.strings["acts.calo.jets.inputEdmHits"] = {m_config.inputEdmHits};
  result.strings["acts.calo.jets.native"] = {m_config.inputJets};
  result.strings["acts.calo.jets.output"] = {m_config.outputJets};
  result.numbers["acts.calo.jets.radius"] = {m_config.reconstruction.radius};
  result.numbers["acts.calo.jets.ptMin"] = {m_config.reconstruction.jetPtMin /
                                            Acts::UnitConstants::GeV};
  const auto& origin = m_config.reconstruction.origin;
  result.numbers["acts.calo.jets.origin"] = {
      origin.x() / Acts::UnitConstants::mm,
      origin.y() / Acts::UnitConstants::mm,
      origin.z() / Acts::UnitConstants::mm};
  return result;
}
}  // namespace ActsExamples
