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

/// @brief This method labels events as either soft or hard and sets the
/// multiplicity
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
      // Search for the maximum in particles with the same PDG ID as the
      // interacting one
      if (p.pdg() == event.initialParticle.pdg())
        maxMom = std::max(p.absMomentum(), maxMom);
      // And the maximum among the others
      else
        maxMomOthers = std::max(p.absMomentum(), maxMomOthers);
    }
    // Label the indication that the interacting particle carries most of the
    // momentum
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

/// @brief This method builds components for transforming a probability distribution into components that represent the cumulative probability distribution
///
/// @param [in] hist The probability distribution
///
/// @return Tuple containing the bin borders, the non-normalised cumulative distribution and the sum over all entries
std::tuple<std::vector<float>, std::vector<double>, double>
buildNotNormalisedMap(TH1F const* hist) {
  // Retrieve the number of bins & borders
  const int nBins = hist->GetNbinsX();
  std::vector<float> histoBorders(nBins + 1);

  // Fill the cumulative histogram
  double integral = 0.;
  std::vector<double> temp_HistoContents(nBins);
  int iBin;
  for (iBin = 0; iBin < nBins; iBin++) {
    float binval = hist->GetBinContent(iBin + 1);
    // Avoid negative bin values
    if (binval < 0) {
      binval = 0.;
    }
    // Store the value
    integral += binval;
    temp_HistoContents[iBin] = integral;
  }

  // Ensure that content is available
  if (integral == 0.) {
    histoBorders.clear();
    temp_HistoContents.clear();
    return std::make_tuple(histoBorders, temp_HistoContents, integral);
  }

  // Set the bin borders
  for (iBin = 1; iBin <= nBins; iBin++)
    histoBorders[iBin - 1] = hist->GetXaxis()->GetBinLowEdge(iBin);
  histoBorders[nBins] = hist->GetXaxis()->GetXmax();

  return std::make_tuple(histoBorders, temp_HistoContents, integral);
}

/// @brief This function combines neighbouring bins with the same value
///
/// @param [in, out] histoBorders The borders of the bins
/// @param [in, out] histoContents The content of each bin
void reduceMap(std::vector<float>& histoBorders, std::vector<uint32_t>& histoContents)
{
	for(auto cit = histoContents.cbegin(); cit != histoContents.cend(); cit++)
	{
		while(std::next(cit, 1) != histoContents.end() && *cit == *std::next(cit, 1))
		{
			const auto distance = std::distance(histoContents.cbegin(), std::next(cit, 1));
			// Remove the bin
			histoBorders.erase(histoBorders.begin() + distance);
			histoContents.erase(histoContents.begin() + distance);
		}
	}
}

/// @brief This method transforms a probability distribution into components
/// that represent the cumulative probability distribution
///
/// @note This method is used to avoid Root dependencies afterwards by
/// decomposing the histogram
/// @param [in] hist The probability distribution
///
/// @return Pair containing the bin borders and the bin content
std::pair<std::vector<float>, std::vector<uint32_t>> buildMap(TH1F const* hist) {
  // Build the components
  std::tuple<std::vector<float>, std::vector<double>, double> map = buildNotNormalisedMap(hist);
  const int nBins = hist->GetNbinsX();
  std::vector<double>& histoContents = std::get<1>(map);
  
  // Fast exit if the histogram is empty
  if(histoContents.empty())
	return std::make_pair(std::get<0>(map), std::vector<uint32_t>());

  // Set the bin content
  std::vector<uint32_t> normalisedHistoContents(nBins);
  const double invIntegral = 1. / std::get<2>(map);
  for (int iBin = 0; iBin < nBins; ++iBin) {
    normalisedHistoContents[iBin] = UINT32_MAX * (histoContents[iBin] * invIntegral);
  }
  
  auto histoBorders = std::get<0>(map);
  reduceMap(histoBorders, normalisedHistoContents);
  return std::make_pair(histoBorders, normalisedHistoContents);
}

/// @brief This method transforms a probability distribution into components
/// that represent the cumulative probability distribution
///
/// @note This method is used to avoid Root dependencies afterwards by
/// decomposing the histogram
/// @param [in] hist The probability distribution
/// @param [in] integral Scaling factor of the distribution
///
/// @return Pair containing the bin borders and the bin content
//~ std::pair<std::vector<float>, std::vector<uint32_t>> buildMap(TH1F const* hist, double integral) {
  //~ // Build the components
  //~ std::tuple<std::vector<float>, std::vector<double>, double> map = buildNotNormalisedMap(hist);

  //~ const int nBins = hist->GetNbinsX();
  //~ std::vector<double>& histoContents = std::get<1>(map);

  //~ // Fast exit if the histogram is empty
  //~ if(histoContents.empty())
	//~ return std::make_pair(std::get<0>(map), std::vector<uint32_t>());

  //~ // Set the bin content
  //~ std::vector<uint32_t> normalisedHistoContents(nBins);
  //~ const double invIntegral = 1. / integral;
  //~ for (int iBin = 0; iBin < nBins; ++iBin) {
    //~ normalisedHistoContents[iBin] = UINT32_MAX * (histoContents[iBin] * invIntegral);
  //~ }

  //~ std::vector<float> histoBorders = std::get<0>(map);
  //~ reduceMap(histoBorders, normalisedHistoContents);

  //~ return std::make_pair(histoBorders, normalisedHistoContents);
//~ }

/// @brief This method builds decomposed cumulative probability distributions
/// out of a vector of proability distributions
///
/// @param [in] histos Vector of probability distributions
///
/// @return Vector containing the decomposed cumulative probability
/// distributions
std::vector<std::pair<std::vector<float>, std::vector<uint32_t>>> buildMaps(
    const std::vector<TH1F*>& histos) {
  std::vector<std::pair<std::vector<float>, std::vector<uint32_t>>> maps;
  for (auto& h : histos) {
    maps.push_back(buildMap(h));
  }
  return maps;
}

/// @brief This method parametrises and stores recursively a parametrisation of
/// the final state kinematics in a nuclear interaction
///
/// @param [in] eventFractionCollection The event storage
/// @param [in] interactionType The interaction type that will be parametrised
/// @param [in] cfg Configuration that steers the binning of histograms
void recordKinematicParametrisation(
    const std::vector<NuclearInteractionParametrisation::EventFraction>&
        eventFractionCollection,
    bool interactionType, unsigned int multiplicity,
    const ActsExamples::RootNuclearInteractionParametersWriter::Config& cfg) {
  gDirectory->mkdir(std::to_string(multiplicity).c_str());
  gDirectory->cd(std::to_string(multiplicity).c_str());

  // Parametrise the momentum und invarian mass distributions
  const auto momentumParameters =
      NuclearInteractionParametrisation::buildMomentumParameters(eventFractionCollection, multiplicity, interactionType,
                          cfg.momentumBins);
  std::vector<NuclearInteractionParametrisation::CumulativeDistribution>
      distributionsMom = momentumParameters.second;
  const auto invariantMassParameters =
      NuclearInteractionParametrisation::buildInvariantMassParameters(eventFractionCollection, multiplicity, interactionType,
                          cfg.invariantMassBins);
  std::vector<NuclearInteractionParametrisation::CumulativeDistribution>
      distributionsInvMass = invariantMassParameters.second;

  // Fast exit in case of no events
  if (!distributionsMom.empty() && !distributionsInvMass.empty()) {
    // Write the eigenspace components for the momenta
    NuclearInteractionParametrisation::EigenspaceComponents
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
    NuclearInteractionParametrisation::EigenspaceComponents
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

	const auto momDistributions = buildMaps(distributionsMom);
	const auto invMassDistributions = buildMaps(distributionsInvMass);

	// Write the distributions
	for (unsigned int i = 0; i <= multiplicity; i++) {
		if(cfg.writeHistograms)
		{
			gDirectory->WriteObject(
			  distributionsMom[i],
			  ("MomentumDistributionHistogram_" + std::to_string(i)).c_str());
	  }
	  gDirectory->WriteObject(&momDistributions[i].first, ("MomentumDistributionBinBorders_" + std::to_string(i)).c_str());
	  gDirectory->WriteObject(&momDistributions[i].second, ("MomentumDistributionBinContents_" + std::to_string(i)).c_str());
	  delete (distributionsMom[i]);
	}
	for (unsigned int i = 0; i < multiplicity; i++) {
		if(cfg.writeHistograms)
		{
		  gDirectory->WriteObject(
			  distributionsInvMass[i],
			  ("InvariantMassDistributionHistogram_" + std::to_string(i)).c_str());
	  }
	  gDirectory->WriteObject(&invMassDistributions[i].first, ("InvariantMassDistributionBinBorders_" + std::to_string(i)).c_str());
	  gDirectory->WriteObject(&invMassDistributions[i].second, ("InvariantMassDistributionBinContents_" + std::to_string(i)).c_str());
	  delete (distributionsInvMass[i]);
	}
  }

  gDirectory->cd("..");
}
}  // namespace

ActsExamples::RootNuclearInteractionParametersWriter::
    RootNuclearInteractionParametersWriter(
        const ActsExamples::RootNuclearInteractionParametersWriter::Config& cfg,
        Acts::Logging::Level lvl)
    : WriterT(cfg.inputSimulationProcesses, "RootNuclearInteractionParametersWriter",
              lvl),
      m_cfg(cfg) {
  if (m_cfg.inputSimulationProcesses.empty()) {
    throw std::invalid_argument("Missing input collection");
  }
  if (m_cfg.outputFilename.empty()) {
    throw std::invalid_argument("Missing output filename for parameters");
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
  if(m_cfg.writeHistograms)
	gDirectory->WriteObject(nuclearInteractionProbability, "NuclearInteractionHistogram");
const auto mapNIprob = buildMap(nuclearInteractionProbability);
gDirectory->WriteObject(&mapNIprob.first, "NuclearInteractionBinBorders");
gDirectory->WriteObject(&mapNIprob.second, "NuclearInteractionBinContents");
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
          m_eventFractionCollection, m_cfg.multiplicityMax);
  gDirectory->cd("soft");
  if(m_cfg.writeHistograms)
	gDirectory->WriteObject(multiplicity.first, "MultiplicityHistogram");
const auto multProbSoft = buildMap(multiplicity.first);
gDirectory->WriteObject(&multProbSoft.first, "MultiplicityBinBorders");
gDirectory->WriteObject(&multProbSoft.second, "MultiplicityBinContents");

  for(unsigned int i = 1; i <= m_cfg.multiplicityMax; i++)
  {
	  recordKinematicParametrisation(m_eventFractionCollection, true, 5, m_cfg);
  }
  gDirectory->cd("../hard");
  if(m_cfg.writeHistograms)
	gDirectory->WriteObject(multiplicity.second, "MultiplicityHistogram");
const auto multProbHard = buildMap(multiplicity.second);
gDirectory->WriteObject(&multProbHard.first, "MultiplicityBinBorders");
gDirectory->WriteObject(&multProbHard.second, "MultiplicityBinContents");

  for(unsigned int i = 1; i <= m_cfg.multiplicityMax; i++)
  {
	recordKinematicParametrisation(m_eventFractionCollection, false, 5, m_cfg);
  }
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
    const std::vector<
        std::tuple<ActsExamples::SimParticle, ActsExamples::SimParticle,
                   std::vector<ActsExamples::SimParticle>>>& event) {
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
