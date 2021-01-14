// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Fatras/NuclearInteractionOptions.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <string>
#include <TFile.h>
#include <TH1F.h>
#include <TVectorF.h>

namespace {
// TODO: some functions were moved to the writer and can be removed
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
std::pair<std::vector<float>, std::vector<uint32_t>> buildMap(TH1F const* hist, double integral) {
  // Build the components
  std::tuple<std::vector<float>, std::vector<double>, double> map = buildNotNormalisedMap(hist);

  const int nBins = hist->GetNbinsX();
  std::vector<double>& histoContents = std::get<1>(map);

  // Fast exit if the histogram is empty
  if(histoContents.empty())
	return std::make_pair(std::get<0>(map), std::vector<uint32_t>());

  // Set the bin content
  std::vector<uint32_t> normalisedHistoContents(nBins);
  const double invIntegral = 1. / integral;
  for (int iBin = 0; iBin < nBins; ++iBin) {
    normalisedHistoContents[iBin] = UINT32_MAX * (histoContents[iBin] * invIntegral);
  }

  std::vector<float> histoBorders = std::get<0>(map);
  reduceMap(histoBorders, normalisedHistoContents);

  return std::make_pair(histoBorders, normalisedHistoContents);}
  
std::vector<std::pair<std::vector<float>, std::vector<uint32_t>>> buildMaps(
    const std::vector<TH1F*>& histos) {
  std::vector<std::pair<std::vector<float>, std::vector<uint32_t>>> maps;
  for (auto& h : histos) {
    maps.push_back(buildMap(h));
  }
  return maps;
}
}

void ActsExamples::Options::addNuclearInteractionOptions(
    ActsExamples::Options::Description& desc) {
  using boost::program_options::value;

  auto opt = desc.add_options();
  opt("fatras-nuclear-interaction-parametrisation",
     value<std::string>()->default_value({}),
     "List of files containing a parametrisation for nuclear interaction.");
  opt("fatras-simulated-events-nuclear-interaction-parametrisation",
    value<int32_t>()->default_value(0),
    "Number of events simulated for the parametrisation of the nuclear interaction.");
}

ActsFatras::detail::MultiParticleParametrisation ActsExamples::Options::readParametrisations(
    const std::string& fileName, const uint32_t& nSimulatedEvents) {
std::cout << "Read call " << std::endl;	
	
	// The collection
	ActsFatras::detail::MultiParticleParametrisation mpp;

    // Now read file
std::cout << "Reading file " << fileName << std::endl;
		ActsFatras::detail::Parametrisation parametrisation; // TODO: Does a single file break the reading pattern?
		TFile tf(fileName.c_str(), "read");
		gDirectory->cd();
		// Walk over all elements in the file
		auto list = gDirectory->GetListOfKeys();
		auto elem = list->First();
		while(elem)
		{
			// Parameters for a fixed inital momentum
			ActsFatras::detail::Parameters parameters;
			// Get the initial momentum
			char const* name = elem->GetName();
			parameters.momentum = std::stof(name);
			gDirectory->cd(name);
std::cout << "Momentum found: "<< parameters.momentum << " " << name << std::endl;
			// Get the nuclear interaction probability
			TH1F* nuclearInteraction = (TH1F*) gDirectory->Get("NuclearInteraction");
std::cout << "NuclearInteraction retrieved: " << nuclearInteraction << std::endl;
			parameters.nuclearInteractionProbability = buildMap(nuclearInteraction, nSimulatedEvents);
std::cout << "Nuclear interaction" << std::endl;
			// Get the soft interaction probability
			TVectorF* softInteraction = (TVectorF*) gDirectory->Get("SoftInteraction");
			parameters.softInteractionProbability = (*softInteraction)[0];
std::cout << "Soft Interaction" << std::endl;
			// Get the branching probabilities
			std::vector<int> branchingPdgIds = *((std::vector<int>*) gDirectory->Get("BranchingPdgIds"));
			std::vector<int> targetPdgIds = *((std::vector<int>*) gDirectory->Get("TargetPdgIds"));
			std::vector<float> targetPdgProbability = *((std::vector<float>*) gDirectory->Get("TargetPdgProbability"));
			parameters.pdgMap.reserve(branchingPdgIds.size());
			for(unsigned int i = 0; i < branchingPdgIds.size(); i++)
			{
				auto it = parameters.pdgMap.begin();
				while(it->first != branchingPdgIds[i] && it != parameters.pdgMap.end())
					it++;
				
				const auto target = std::make_pair(targetPdgIds[i], targetPdgProbability[i]);
				if(it != parameters.pdgMap.end())
				{
					it->second.push_back(target);
				} else {
					parameters.pdgMap.push_back(std::make_pair(branchingPdgIds[i], std::vector<std::pair<int, float>>{target}));
				}
			}
std::cout << "pdg probability" << std::endl;
			// Get the soft distributions
			gDirectory->cd("soft");
			// Get the multiplicity distribution
			TH1F* softMultiplicity = (TH1F*) gDirectory->Get("Multiplicity");
			parameters.softMultiplicity = buildMap(softMultiplicity);
std::cout << "multiplicity" << std::endl;
			//~ // Get the distributions for each final state multiplicity
			auto softList = gDirectory->GetListOfKeys();
			auto softElement = softList->First();
			while(softElement)
			{
				if(softElement->IsFolder())
				{
					// Find the momentum and invariant mass distributions
					const char* distributionName = softElement->GetName();
					unsigned int mult = std::stoi(distributionName);
					gDirectory->cd(distributionName);
					std::vector<TH1F*> momentumDistributions;
					momentumDistributions.reserve(mult + 1);
					std::vector<TH1F*> invariantMassDistributions;
					invariantMassDistributions.reserve(mult);
					for(unsigned int i = 0; i < mult; i++)
					{
						momentumDistributions.push_back(std::move((TH1F*) gDirectory->Get(("MomentumDistribution_" + std::to_string(i)).c_str())));
						invariantMassDistributions.push_back(std::move((TH1F*) gDirectory->Get(("InvariantMassDistribution_" + std::to_string(i)).c_str())));
					}
					momentumDistributions.push_back(std::move((TH1F*) gDirectory->Get(("MomentumDistribution_" + std::to_string(mult)).c_str())));

					// Get the eigenspace components for the kinematic parameters
					std::vector<float> momentumEigenvalues = *((std::vector<float>*) gDirectory->Get("MomentumEigenvalues"));
					std::vector<float> momentumEigenvectors = *((std::vector<float>*) gDirectory->Get("MomentumEigenvectors"));
					std::vector<float> momentumMean = *((std::vector<float>*) gDirectory->Get("MomentumMean"));
					std::vector<float> invariantMassEigenvalues = *((std::vector<float>*) gDirectory->Get("InvariantMassEigenvalues"));
					std::vector<float> invariantMassEigenvectors = *((std::vector<float>*) gDirectory->Get("InvariantMassEigenvectors"));
					std::vector<float> invariantMassMean = *((std::vector<float>*) gDirectory->Get("InvariantMassMean"));
					
					// Prepare the storage
					if(mult >= parameters.softKinematicParameters.size())
						parameters.softKinematicParameters.resize(mult + 1);
					// Prepare and store the kinematic parameters
					ActsFatras::detail::Parameters::ParametersWithFixedMultiplicity kinematicParameters(buildMaps(momentumDistributions), 
					momentumEigenvalues, momentumEigenvectors, momentumMean,
					buildMaps(invariantMassDistributions),
					invariantMassEigenvalues, invariantMassEigenvectors, invariantMassMean);
					parameters.softKinematicParameters[mult] = kinematicParameters;
					
					gDirectory->cd("..");
				}
				softElement = softList->After(softElement);
			}
			// Get the hard distributions
			gDirectory->cd("../hard");
			// Get the multiplicity distribution
			TH1F* hardMultiplicity = (TH1F*) gDirectory->Get("Multiplicity");
			parameters.hardMultiplicity = buildMap(hardMultiplicity);

			// Get the distributions for each final state multiplicity
			auto hardList = gDirectory->GetListOfKeys();
			auto hardElement = hardList->First();
			while(hardElement)
			{
				if(hardElement->IsFolder())
				{
					// Find the momentum and invariant mass distributions
					const char* distributionName = hardElement->GetName();
					unsigned int mult = std::stoi(distributionName);
					gDirectory->cd(distributionName);
					std::vector<TH1F*> momentumDistributions;
					momentumDistributions.reserve(mult + 1);
					std::vector<TH1F*> invariantMassDistributions;
					invariantMassDistributions.reserve(mult);
					for(unsigned int i = 0; i < mult; i++)
					{
						momentumDistributions.push_back(std::move((TH1F*) gDirectory->Get(("MomentumDistribution_" + std::to_string(i)).c_str())));
						invariantMassDistributions.push_back(std::move((TH1F*) gDirectory->Get(("InvariantMassDistribution_" + std::to_string(i)).c_str())));
					}
					momentumDistributions.push_back(std::move((TH1F*) gDirectory->Get(("MomentumDistribution_" + std::to_string(mult)).c_str())));
	
					// Get the eigenspace components for the kinematic parameters
					std::vector<float> momentumEigenvalues = *((std::vector<float>*) gDirectory->Get("MomentumEigenvalues"));
					std::vector<float> momentumEigenvectors = *((std::vector<float>*) gDirectory->Get("MomentumEigenvectors"));
					std::vector<float> momentumMean = *((std::vector<float>*) gDirectory->Get("MomentumMean"));
					std::vector<float> invariantMassEigenvalues = *((std::vector<float>*) gDirectory->Get("InvariantMassEigenvalues"));
					std::vector<float> invariantMassEigenvectors = *((std::vector<float>*) gDirectory->Get("InvariantMassEigenvectors"));
					std::vector<float> invariantMassMean = *((std::vector<float>*) gDirectory->Get("InvariantMassMean"));
						
					// Prepare the storage			
					if(mult >= parameters.hardKinematicParameters.size())
						parameters.hardKinematicParameters.resize(mult + 1);
					// Prepare and store the kinematic parameters
					ActsFatras::detail::Parameters::ParametersWithFixedMultiplicity kinematicParameters(buildMaps(momentumDistributions), 
					momentumEigenvalues, momentumEigenvectors, momentumMean,
					buildMaps(invariantMassDistributions),
					invariantMassEigenvalues, invariantMassEigenvectors, invariantMassMean);
					parameters.hardKinematicParameters[mult] = kinematicParameters;

					gDirectory->cd("..");
				}
				hardElement = softList->After(hardElement);
			}
			elem = list->After(elem); // TODO: this might be not needed
			// Store the parametrisation
			parametrisation.push_back(std::make_pair(parameters.momentum, parameters));
		}
		tf.Close();
		
		// Write to the collection to the EventStore
		mpp.push_back(std::make_pair(211, parametrisation));
  
  // Return success flag
  return mpp;
}
