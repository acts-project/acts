// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootNuclearInteractionParametersWriter.hpp"

#include "ActsExamples/EventData/SimParticle.hpp"
#include <stdexcept>
#include <TFile.h>

namespace {

/// @brief This method labels events as either soft or hard and sets the multiplicity
///
/// @param [in] events The events that will be labeled
void labelEvents(
    std::vector<NuclearInteractionParametrisation::EventFraction>& events) {
  // Search for the highest momentum particles per event
  for (NuclearInteractionParametrisation::EventFraction& event : events) {
    double maxMom = 0.;
    double maxMomOthers = 0.;
    // Walk over all final state particles
    for (const ActsExamples::SimParticle& p : event.finalParticles) {
		// Search for the maximum in particles with the same PDG ID as the interacting one
      if (p.pdg() == event.initialParticle.pdg())
        maxMom = std::max(p.absMomentum(), maxMom);
        // And the maximum among the others
      else
        maxMomOthers = std::max(p.absMomentum(), maxMomOthers);
    }
    // Label the indication that the interacting particle carries most of the momentum
    event.soft = (maxMom > maxMomOthers);

	// Get the final state p_T
    double pt = 0.;
    Acts::Vector2D ptVec(0., 0.);
    for (const ActsExamples::SimParticle& p : event.finalParticles) {
      Acts::Vector2D particlePt =
          p.momentum4().template segment<2>(Acts::eMom0);
      ptVec[0] += particlePt[0];
      ptVec[1] += particlePt[1];
    }
    pt = ptVec.norm();

	// Use the final state p_T as veto for the soft label
    if (event.soft && pt <= event.interactingParticle.transverseMomentum())
      event.soft = false;

	// Store the multiplicity
    event.multiplicity = event.finalParticles.size();
  }
}

/// @brief This method parametrises and stores recursively a parametrisation of
/// the final state kinematics in a nuclear interaction
///
/// @tparam multiplicity_t The final state multiplicity
/// @param [in] eventFractionCollection The event storage
/// @param [in] interactionType The interaction type that will be parametrised
/// @param [in] cfg Configuration that steers the binning of histograms
template <unsigned int multiplicity_t>
void recordKinematicParametrisation(
    const std::vector<NuclearInteractionParametrisation::EventFraction>&
        eventFractionCollection,
    bool interactionType,
    const ActsExamples::RootNuclearInteractionParametersWriter::Config& cfg) {
  gDirectory->mkdir(std::to_string(multiplicity_t).c_str());
  gDirectory->cd(std::to_string(multiplicity_t).c_str());

  // Parametrise the momentum und invarian mass distributions
  const auto momentumParameters =
      NuclearInteractionParametrisation::buildMomentumParameters<
          multiplicity_t>(eventFractionCollection, interactionType,
                          cfg.momentumBins);
  std::vector<NuclearInteractionParametrisation::CumulativeDistribution>
      distributionsMom = momentumParameters.second;
  const auto invariantMassParameters =
      NuclearInteractionParametrisation::buildInvariantMassParameters<
          multiplicity_t>(eventFractionCollection, interactionType,
                          cfg.invariantMassBins);
  std::vector<NuclearInteractionParametrisation::CumulativeDistribution>
      distributionsInvMass = invariantMassParameters.second;

  // Fast exit in case of no events
  if (!distributionsMom.empty() && !distributionsInvMass.empty()) {
    // Write the eigenspace components for the momenta
    NuclearInteractionParametrisation::EigenspaceComponents<multiplicity_t + 1>
        esComponentsMom = momentumParameters.first;

    auto momEigenVal = std::get<0>(esComponentsMom);
    auto momEigenVec = std::get<1>(esComponentsMom);
    auto momMean = std::get<2>(esComponentsMom);
    std::vector<float> momVecVal(momEigenVal.data(),
                                 momEigenVal.data() + momEigenVal.size());
    std::vector<float> momVecVec(momEigenVec.data(),
                                 momEigenVec.data() + momEigenVec.size());
    std::vector<float> momVecMean(momMean.data(),
                                  momMean.data() + momMean.size());

    gDirectory->WriteObject(&momVecVal, "MomentumEigenvalues");
    gDirectory->WriteObject(&momVecVec, "MomentumEigenvectors");
    gDirectory->WriteObject(&momVecMean, "MomentumMean");

    // Write the eigenspace components for the invariant masses
    NuclearInteractionParametrisation::EigenspaceComponents<multiplicity_t + 1>
        esComponentsInvMass = invariantMassParameters.first;

    auto invMassEigenVal = std::get<0>(esComponentsInvMass);
    auto invMassEigenVec = std::get<1>(esComponentsInvMass);
    auto invMassMean = std::get<2>(esComponentsInvMass);
    std::vector<float> invMassVecVal(
        invMassEigenVal.data(),
        invMassEigenVal.data() + invMassEigenVal.size());
    std::vector<float> invMassVecVec(
        invMassEigenVec.data(),
        invMassEigenVec.data() + invMassEigenVec.size());
    std::vector<float> invMassVecMean(invMassMean.data(),
                                      invMassMean.data() + invMassMean.size());

    gDirectory->WriteObject(&invMassVecVal, "InvariantMassEigenvalues");
    gDirectory->WriteObject(&invMassVecVec, "InvariantMassEigenvectors");
    gDirectory->WriteObject(&invMassVecMean, "InvariantMassMean");

    // Write the distributions
    for (unsigned int i = 0; i <= multiplicity_t; i++) {
      gDirectory->WriteObject(
          distributionsMom[i],
          ("MomentumDistribution_" + std::to_string(i)).c_str());
      delete (distributionsMom[i]);
    }
    for (unsigned int i = 0; i < multiplicity_t; i++) {
      gDirectory->WriteObject(
          distributionsInvMass[i],
          ("InvariantMassDistribution_" + std::to_string(i)).c_str());
      delete (distributionsInvMass[i]);
    }
  }

  gDirectory->cd("..");
  // Repeat for (multiplicity - 1)
  recordKinematicParametrisation<multiplicity_t - 1>(eventFractionCollection,
                                                     interactionType, cfg);
}
/// Recursion break
template <>
void recordKinematicParametrisation<0>(
    const std::vector<NuclearInteractionParametrisation::
                          EventFraction>& /*eventFractionCollection*/,
    bool /*interactionType*/,
    const ActsExamples::RootNuclearInteractionParametersWriter::
        Config& /*cfg*/) {}
}  // namespace

ActsExamples::RootNuclearInteractionParametersWriter::
    RootNuclearInteractionParametersWriter(
        const ActsExamples::RootNuclearInteractionParametersWriter::Config& cfg,
        Acts::Logging::Level lvl)
    : WriterT(cfg.inputEventFractions, "RootNuclearInteractionParametersWriter",
              lvl),
      m_cfg(cfg) {
  if (m_cfg.inputEventFractions.empty()) {
    throw std::invalid_argument("Missing input collection");
  }
  if (m_cfg.outputFilename.empty()) {
    throw std::invalid_argument("Missing output filename");
  }
}

ActsExamples::RootNuclearInteractionParametersWriter::
    ~RootNuclearInteractionParametersWriter() {}

ActsExamples::ProcessCode
ActsExamples::RootNuclearInteractionParametersWriter::endRun() {
  if (m_eventFractionCollection.empty())
    return ProcessCode::ABORT;

  // Exclusive access to the file while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // The file
  TFile tf(m_cfg.outputFilename.c_str(), m_cfg.fileMode.c_str());
  gDirectory->cd();
  gDirectory->mkdir(
      std::to_string(m_eventFractionCollection[0].initialParticle.absMomentum())
          .c_str());
  gDirectory->cd(
      std::to_string(m_eventFractionCollection[0].initialParticle.absMomentum())
          .c_str());
  gDirectory->mkdir("soft");
  gDirectory->mkdir("hard");

  // Write the nuclear interaction probability
  const auto nuclearInteractionProbability = NuclearInteractionParametrisation::
      cumulativeNuclearInteractionProbability(m_eventFractionCollection,
                                              m_cfg.interactionProbabilityBins);
  gDirectory->WriteObject(nuclearInteractionProbability, "NuclearInteraction");
  delete (nuclearInteractionProbability);

  // Write the interaction type proability
  const auto softProbability =
      NuclearInteractionParametrisation::softProbability(
          m_eventFractionCollection);
  gDirectory->WriteObject(&softProbability, "SoftInteraction");

  // Write the PDG id production distribution
  const auto pdgIdMap =
      NuclearInteractionParametrisation::cumulativePDGprobability(
          m_eventFractionCollection);
  std::vector<int> branchingPdgIds;
  std::vector<int> targetPdgIds;
  std::vector<float> targetPdgProbability;
  for (const auto& targetPdgIdMap : pdgIdMap) {
    for (const auto& producedPdgIdMap : targetPdgIdMap.second) {
      branchingPdgIds.push_back(targetPdgIdMap.first);
      targetPdgIds.push_back(producedPdgIdMap.first);
      targetPdgProbability.push_back(producedPdgIdMap.second);
    }
  }
  gDirectory->WriteObject(&branchingPdgIds, "BranchingPdgIds");
  gDirectory->WriteObject(&targetPdgIds, "TargetPdgIds");
  gDirectory->WriteObject(&targetPdgProbability, "TargetPdgProbability");

  // Write the multiplicity and kinematics distribution
  const auto multiplicity =
      NuclearInteractionParametrisation::cumulativeMultiplicityProbability(
          m_eventFractionCollection);
  gDirectory->cd("soft");
  gDirectory->WriteObject(multiplicity.first, "Multiplicity");
  recordKinematicParametrisation<5>(m_eventFractionCollection, true, m_cfg);
  gDirectory->cd("../hard");
  gDirectory->WriteObject(multiplicity.second, "Multiplicity");
  recordKinematicParametrisation<5>(m_eventFractionCollection, false, m_cfg);
  delete (multiplicity.first);
  delete (multiplicity.second);

  gDirectory->cd();
  tf.Write();
  tf.Close();

  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode
ActsExamples::RootNuclearInteractionParametersWriter::writeT(
    const AlgorithmContext& /*ctx*/,
    const std::vector<std::tuple<ActsExamples::SimParticle, ActsExamples::SimParticle,
                                std::vector<ActsExamples::SimParticle>>>&
        event) {
  // Convert the tuple to use additional categorisation variables
  std::vector<NuclearInteractionParametrisation::EventFraction> eventFractions;
  eventFractions.reserve(event.size());
  for (const auto& e : event) {
    eventFractions.emplace_back(e);
  }
  // Categorise the event
  labelEvents(eventFractions);
  // Store the event
  m_eventFractionCollection.insert(m_eventFractionCollection.end(),
                                   eventFractions.begin(),
                                   eventFractions.end());

  return ProcessCode::SUCCESS;
}
