// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Calorimeter/CalorimeterClusteringAlgorithm.hpp"
#include "ActsExamples/Calorimeter/CalorimeterDigitizationAlgorithm.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"

#ifdef ACTS_EXAMPLE_CALORIMETER_JETS
#include "ActsExamples/Calorimeter/CalorimeterJetAlgorithm.hpp"
#endif

#include <cmath>
#include <memory>
#include <utility>

namespace {

using namespace Acts::UnitLiterals;

class SyntheticDeposits final : public ActsExamples::IAlgorithm {
 public:
  SyntheticDeposits()
      : IAlgorithm("SyntheticCaloDeposits",
                   Acts::getDefaultLogger("SyntheticCaloDeposits",
                                          Acts::Logging::INFO)) {
    m_output.initialize("sim_calo_hits");
  }

  ActsExamples::ProcessCode execute(
      const ActsExamples::AlgorithmContext& context) const override {
    const Acts::Vector3 centre{1500_mm, 0, 0};
    // Two deposits merge into cell 42, cell 43 neighbours it, cell 46 is
    // isolated. Cell 44 is below threshold and cell 45 is outside the time
    // window. No detector or shower is simulated.
    ActsCalorimeter::SimCalorimeterHitContainer deposits{
        {42, centre, 0.01_GeV, 1_ns},
        {42, centre, 0.02_GeV, 4_ns},
        {43, Acts::Vector3{1500_mm, 10_mm, 0}, 0.02_GeV, 5_ns},
        {44, centre, 0.001_GeV, 2_ns},
        {45, centre, 1_GeV, 100_ns},
        {46, Acts::Vector3{-1500_mm, 0, 0}, 0.03_GeV, 2_ns}};
    m_output(context, std::move(deposits));
    return ActsExamples::ProcessCode::SUCCESS;
  }

 private:
  ActsExamples::WriteDataHandle<ActsCalorimeter::SimCalorimeterHitContainer>
      m_output{this, "OutputSimHits"};
};

class CheckCalibratedHits final : public ActsExamples::IAlgorithm {
 public:
  CheckCalibratedHits()
      : IAlgorithm("CheckCaloHits", Acts::getDefaultLogger(
                                        "CheckCaloHits", Acts::Logging::INFO)) {
    m_input.initialize("calo_hits");
  }

  ActsExamples::ProcessCode execute(
      const ActsExamples::AlgorithmContext& context) const override {
    const auto& hits = m_input(context);
    if (hits.size() != 3 || hits.front().cellId != 42 ||
        std::abs(hits.front().energy - 0.06_GeV) > 1e-12_GeV ||
        std::abs(hits.front().time - 3_ns) > 1e-12_ns ||
        hits.front().sourceIndices != std::vector<std::size_t>{0, 1} ||
        hits[1].cellId != 43 ||
        std::abs(hits[1].energy - 0.04_GeV) > 1e-12_GeV ||
        hits[2].cellId != 46 ||
        std::abs(hits[2].energy - 0.06_GeV) > 1e-12_GeV) {
      ACTS_ERROR("Unexpected prototype calorimeter response");
      return ActsExamples::ProcessCode::ABORT;
    }
    ACTS_INFO("Event " << context.eventNumber << ": cell 42, energy "
                       << hits.front().energy / 1_GeV << " GeV, time "
                       << hits.front().time / 1_ns << " ns");
    return ActsExamples::ProcessCode::SUCCESS;
  }

 private:
  ActsExamples::ReadDataHandle<ActsCalorimeter::CalorimeterHitContainer>
      m_input{this, "InputHits"};
};

class CheckClusters final : public ActsExamples::IAlgorithm {
 public:
  CheckClusters()
      : IAlgorithm(
            "CheckCaloClusters",
            Acts::getDefaultLogger("CheckCaloClusters", Acts::Logging::INFO)) {
    m_input.initialize("calo_clusters");
  }

  ActsExamples::ProcessCode execute(
      const ActsExamples::AlgorithmContext& context) const override {
    const auto& clusters = m_input(context);
    if (clusters.size() != 2 || clusters[0].seedCellId != 42 ||
        std::abs(clusters[0].energy - 0.1_GeV) > 1e-12_GeV ||
        std::abs(clusters[0].position.y() - 4_mm) > 1e-12_mm ||
        std::abs(clusters[0].time - 3.8_ns) > 1e-12_ns ||
        clusters[0].hitIndices != std::vector<std::size_t>{0, 1} ||
        clusters[1].seedCellId != 46 ||
        std::abs(clusters[1].energy - 0.06_GeV) > 1e-12_GeV ||
        clusters[1].hitIndices != std::vector<std::size_t>{2}) {
      ACTS_ERROR("Unexpected prototype calorimeter clusters");
      return ActsExamples::ProcessCode::ABORT;
    }
    ACTS_INFO("Event " << context.eventNumber << ": two clusters with total "
                       << (clusters[0].energy + clusters[1].energy) / 1_GeV
                       << " GeV");
    return ActsExamples::ProcessCode::SUCCESS;
  }

 private:
  ActsExamples::ReadDataHandle<ActsCalorimeter::CalorimeterClusterContainer>
      m_input{this, "InputClusters"};
};

#ifdef ACTS_EXAMPLE_CALORIMETER_JETS
class CheckJets final : public ActsExamples::IAlgorithm {
 public:
  CheckJets()
      : IAlgorithm("CheckCaloJets", Acts::getDefaultLogger(
                                        "CheckCaloJets", Acts::Logging::INFO)) {
    m_input.initialize("calo_jets");
  }

  ActsExamples::ProcessCode execute(
      const ActsExamples::AlgorithmContext& context) const override {
    const auto& jets = m_input(context);
    // The two clusters face opposite directions, so anti-kt R=0.4 keeps them
    // in separate jets. All clusters are transverse to the beam in this toy.
    const double phi = std::atan2(4_mm, 1500_mm);
    const Acts::Vector4 leading{0.1_GeV * std::cos(phi),
                                0.1_GeV * std::sin(phi), 0, 0.1_GeV};
    const Acts::Vector4 subleading{-0.06_GeV, 0, 0, 0.06_GeV};
    if (jets.size() != 2 ||
        (jets[0].fourMomentum - leading).norm() > 1e-12_GeV ||
        (jets[1].fourMomentum - subleading).norm() > 1e-12_GeV ||
        jets[0].clusterIndices != std::vector<std::size_t>{0} ||
        jets[1].clusterIndices != std::vector<std::size_t>{1}) {
      ACTS_ERROR("Unexpected prototype calorimeter jets");
      return ActsExamples::ProcessCode::ABORT;
    }
    ACTS_INFO("Event " << context.eventNumber << ": two calorimeter jets, pT "
                       << jets[0].fourMomentum.head<2>().norm() / 1_GeV << ", "
                       << jets[1].fourMomentum.head<2>().norm() / 1_GeV
                       << " GeV");
    return ActsExamples::ProcessCode::SUCCESS;
  }

 private:
  ActsExamples::ReadDataHandle<ActsCalorimeter::CalorimeterJetContainer>
      m_input{this, "InputJets"};
};
#endif

}  // namespace

int main() {
  ActsExamples::Sequencer::Config sequence;
  sequence.events = 3;
  sequence.numThreads = 1;
  ActsExamples::Sequencer sequencer(sequence);
  sequencer.addAlgorithm(std::make_shared<SyntheticDeposits>());

  ActsExamples::CalorimeterDigitizationAlgorithm::Config digitization;
  digitization.inputSimHits = "sim_calo_hits";
  digitization.outputHits = "calo_hits";
  digitization.response.energyScale = 2;
  digitization.response.energyThreshold = 0.01_GeV;
  digitization.response.timeMin = 0_ns;
  digitization.response.timeMax = 10_ns;
  sequencer.addAlgorithm(
      std::make_shared<ActsExamples::CalorimeterDigitizationAlgorithm>(
          digitization));
  sequencer.addAlgorithm(std::make_shared<CheckCalibratedHits>());
  ActsExamples::CalorimeterClusteringAlgorithm::Config clustering;
  clustering.inputHits = "calo_hits";
  clustering.outputClusters = "calo_clusters";
  clustering.clustering.seedEnergyThreshold = 0.05_GeV;
  clustering.clustering.neighbourEnergyThreshold = 0.02_GeV;
  clustering.neighbours = {{42, 43}};
  sequencer.addAlgorithm(
      std::make_shared<ActsExamples::CalorimeterClusteringAlgorithm>(
          clustering));
  sequencer.addAlgorithm(std::make_shared<CheckClusters>());
#ifdef ACTS_EXAMPLE_CALORIMETER_JETS
  ActsExamples::CalorimeterJetAlgorithm::Config jets;
  jets.inputClusters = "calo_clusters";
  jets.outputJets = "calo_jets";
  jets.jets.radius = 0.4;
  jets.jets.jetPtMin = 0_GeV;
  sequencer.addAlgorithm(
      std::make_shared<ActsExamples::CalorimeterJetAlgorithm>(jets));
  sequencer.addAlgorithm(std::make_shared<CheckJets>());
#endif
  return sequencer.run();
}
