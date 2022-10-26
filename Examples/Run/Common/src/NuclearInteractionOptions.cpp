// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Options/NuclearInteractionOptions.hpp"

#include "ActsExamples/Utilities/Options.hpp"

#include <string>

#include <TFile.h>
#include <TH1F.h>
#include <TVectorF.h>

namespace {
std::pair<std::vector<float>, std::vector<uint32_t>> readHistogram(
    const std::string&& binBordersName, const std::string&& binContentsName) {
  std::vector<float>* binBorders;
  std::vector<uint32_t>* binContents;
  // Get the decomposed histogram
  binBorders = (std::vector<float>*)gDirectory->Get(binBordersName.c_str());
  binContents =
      (std::vector<uint32_t>*)gDirectory->Get(binContentsName.c_str());
  // Return the histogram if available
  if (binBorders != nullptr && binContents != nullptr) {
    return std::make_pair(*binBorders, *binContents);
  }
  return std::make_pair(std::vector<float>(), std::vector<uint32_t>());
}

Acts::ActsDynamicVector readVector(const std::string&& vectorName) {
  std::vector<float>* vector;
  // Get the vector
  vector = (std::vector<float>*)gDirectory->Get(vectorName.c_str());
  // Return the vector if available
  if (vector != nullptr) {
    Acts::ActsDynamicVector result;
    const unsigned int sizeVec = vector->size();
    result.resize(sizeVec);
    for (unsigned int i = 0; i < sizeVec; i++) {
      result(i) = (*vector)[i];
    }

    return result;
  }
  return Acts::ActsDynamicVector();
}

void readKinematicParameters(
    ActsFatras::detail::NuclearInteractionParameters& parameters,
    TObject* folder, bool softInteractionParameters) {
  if (folder->IsFolder()) {
    // Find the momentum and invariant mass distributions
    const char* distributionName = folder->GetName();
    unsigned int mult = std::stoi(distributionName);
    gDirectory->cd(distributionName);
    std::vector<std::pair<std::vector<float>, std::vector<uint32_t>>>
        momentumDistributions;
    momentumDistributions.resize(mult + 1);
    std::vector<std::pair<std::vector<float>, std::vector<uint32_t>>>
        invariantMassDistributions;
    invariantMassDistributions.resize(mult);
    for (unsigned int i = 0; i < mult; i++) {
      momentumDistributions[i] = readHistogram(
          ("MomentumDistributionBinBorders_" + std::to_string(i)).c_str(),
          ("MomentumDistributionBinContents_" + std::to_string(i)).c_str());

      invariantMassDistributions[i] = readHistogram(
          ("InvariantMassDistributionBinBorders_" + std::to_string(i)).c_str(),
          ("InvariantMassDistributionBinContents_" + std::to_string(i))
              .c_str());
    }
    momentumDistributions.back() = readHistogram(
        ("MomentumDistributionBinBorders_" + std::to_string(mult)).c_str(),
        ("MomentumDistributionBinContents_" + std::to_string(mult)).c_str());

    // Get the eigenspace components for the kinematic parameters
    Acts::ActsDynamicVector momentumEigenvalues =
        readVector("MomentumEigenvalues");
    Acts::ActsDynamicVector momentumEigenvectors =
        readVector("MomentumEigenvectors");
    Acts::ActsDynamicVector momentumMean = readVector("MomentumMean");
    Acts::ActsDynamicVector invariantMassEigenvalues =
        readVector("InvariantMassEigenvalues");
    Acts::ActsDynamicVector invariantMassEigenvectors =
        readVector("InvariantMassEigenvectors");
    Acts::ActsDynamicVector invariantMassMean = readVector("InvariantMassMean");

    // Test that a parametrisation is present
    if (momentumEigenvalues.size() != 0 && momentumEigenvectors.size() != 0 &&
        momentumMean.size() != 0 && invariantMassEigenvalues.size() != 0 &&
        invariantMassEigenvectors.size() != 0 &&
        invariantMassMean.size() != 0) {
      // Prepare and store the kinematic parameters
      ActsFatras::detail::NuclearInteractionParameters::
          ParametersWithFixedMultiplicity kinematicParameters(
              momentumDistributions, momentumEigenvalues, momentumEigenvectors,
              momentumMean, invariantMassDistributions,
              invariantMassEigenvalues, invariantMassEigenvectors,
              invariantMassMean);
      if (softInteractionParameters) {
        if (mult >= parameters.softKinematicParameters.size()) {
          parameters.softKinematicParameters.resize(mult + 1);
        }
        parameters.softKinematicParameters[mult] = kinematicParameters;
      } else {
        if (mult >= parameters.hardKinematicParameters.size()) {
          parameters.hardKinematicParameters.resize(mult + 1);
        }
        parameters.hardKinematicParameters[mult] = kinematicParameters;
      }
    }
    gDirectory->cd("..");
  }
}
}  // namespace

void ActsExamples::Options::addNuclearInteractionOptions(
    ActsExamples::Options::Description& desc) {
  using boost::program_options::value;

  auto opt = desc.add_options();
  opt("fatras-nuclear-interaction-parametrisation",
      value<std::string>()->default_value({}),
      "File containing parametrisations for nuclear interaction.");
}

std::string ActsExamples::Options::readNuclearInteractionConfig(
    const boost::program_options::variables_map& variables) {
  return variables["fatras-nuclear-interaction-parametrisation"]
      .as<std::string>();
}

ActsFatras::detail::MultiParticleNuclearInteractionParametrisation
ActsExamples::Options::readParametrisations(const std::string& fileName) {
  // The collection
  ActsFatras::detail::MultiParticleNuclearInteractionParametrisation mpp;

  // Now read file
  ActsFatras::detail::NuclearInteractionParametrisation parametrisation;
  TFile tf(fileName.c_str(), "read");
  gDirectory->cd();
  auto listOfParticles = gDirectory->GetListOfKeys();
  auto initialParticle = listOfParticles->First();
  while (initialParticle != nullptr) {
    // Get the initial particle
    char const* particleName = initialParticle->GetName();
    gDirectory->cd(particleName);

    // Walk over all initial momenta for a particle
    auto listOfMomenta = gDirectory->GetListOfKeys();
    auto initialMomentum = listOfMomenta->First();
    while (initialMomentum != nullptr) {
      // Parameters for a fixed inital momentum
      ActsFatras::detail::NuclearInteractionParameters parameters;
      // Get the initial momentum
      char const* nameMomentum = initialMomentum->GetName();
      parameters.momentum = std::stof(nameMomentum);
      gDirectory->cd(nameMomentum);

      // Get the nuclear interaction probability
      parameters.nuclearInteractionProbability = readHistogram(
          "NuclearInteractionBinBorders", "NuclearInteractionBinContents");

      // Get the soft interaction probability
      TVectorF* softInteraction = (TVectorF*)gDirectory->Get("SoftInteraction");
      parameters.softInteractionProbability = (*softInteraction)[0];

      // Get the branching probabilities
      std::vector<int> branchingPdgIds =
          *((std::vector<int>*)gDirectory->Get("BranchingPdgIds"));
      std::vector<int> targetPdgIds =
          *((std::vector<int>*)gDirectory->Get("TargetPdgIds"));
      std::vector<float> targetPdgProbability =
          *((std::vector<float>*)gDirectory->Get("TargetPdgProbability"));
      parameters.pdgMap.reserve(branchingPdgIds.size());
      for (unsigned int i = 0; i < branchingPdgIds.size(); i++) {
        auto it = parameters.pdgMap.begin();
        while (it->first != branchingPdgIds[i] &&
               it != parameters.pdgMap.end()) {
          it++;
        }

        const auto target =
            std::make_pair(targetPdgIds[i], targetPdgProbability[i]);
        if (it != parameters.pdgMap.end()) {
          it->second.push_back(target);
        } else {
          parameters.pdgMap.push_back(std::make_pair(
              branchingPdgIds[i], std::vector<std::pair<int, float>>{target}));
        }
      }

      // Get the soft distributions
      gDirectory->cd("soft");
      // Get the multiplicity distribution
      parameters.softMultiplicity =
          readHistogram("MultiplicityBinBorders", "MultiplicityBinContents");
      // Get the distributions for each final state multiplicity
      auto softList = gDirectory->GetListOfKeys();
      auto softElement = softList->First();
      while (softElement != nullptr) {
        readKinematicParameters(parameters, softElement, true);
        softElement = softList->After(softElement);
      }

      // Get the hard distributions
      gDirectory->cd("../hard");
      // Get the multiplicity distribution
      parameters.hardMultiplicity =
          readHistogram("MultiplicityBinBorders", "MultiplicityBinContents");

      // Get the distributions for each final state multiplicity
      auto hardList = gDirectory->GetListOfKeys();
      auto hardElement = hardList->First();
      while (hardElement != nullptr) {
        readKinematicParameters(parameters, hardElement, false);
        hardElement = hardList->After(hardElement);
      }

      initialMomentum = listOfMomenta->After(initialMomentum);
      // Store the parametrisation
      parametrisation.push_back(
          std::make_pair(parameters.momentum, parameters));
    }
    tf.Close();

    // Write to the collection to the EventStore
    mpp.push_back(std::make_pair(std::stof(particleName), parametrisation));

    initialParticle = listOfParticles->After(initialParticle);
  }
  // Return success flag
  return mpp;
}
