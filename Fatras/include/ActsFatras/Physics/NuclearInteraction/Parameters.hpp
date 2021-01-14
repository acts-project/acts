// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Common.hpp"

#include <map>
#include <unordered_map>
#include <vector>

namespace ActsFatras {
namespace detail {

/// @brief Data storage of the parametrized nuclear interaction
struct Parameters {
  using CumulativeDistribution =
      std::pair<std::vector<float>, std::vector<uint32_t>>;
  using Distributions = std::vector<CumulativeDistribution>;
  using PdgMap =
      std::vector<std::pair<int, std::vector<std::pair<int, float>>>>;

  /// @brief Nested struct for the storage of the kinematic parametrisations for
  /// a given final state multiplicity
  struct ParametersWithFixedMultiplicity {
    ParametersWithFixedMultiplicity() = default;

    ParametersWithFixedMultiplicity(
        Distributions&& momenta, std::vector<float>& eValMom,
        std::vector<float>& eVecMom, std::vector<float>& meanMom,
        Distributions&& invariantMasses, std::vector<float>& eValIM,
        std::vector<float>& eVecIM, std::vector<float>& meanIM)
        : momentumDistributions(momenta),
          invariantMassDistributions(invariantMasses) {
      const unsigned int sizeMom = eValMom.size();
      eigenvaluesMomentum.resize(sizeMom);
      eigenvectorsMomentum.resize(sizeMom, sizeMom);
      meanMomentum.resize(sizeMom);

      for (unsigned int i = 0; i < sizeMom; i++) {
        eigenvaluesMomentum(i) = eValMom[i];
        for (unsigned int j = 0; j < sizeMom; j++)
          eigenvectorsMomentum(i, j) = eVecMom[i * sizeMom + j];
        meanMomentum(i) = meanMom[i];
      }

      const unsigned int sizeInvMass = eValIM.size();
      eigenvaluesInvariantMass.resize(sizeInvMass);
      eigenvectorsInvariantMass.resize(sizeInvMass, sizeInvMass);
      meanInvariantMass.resize(sizeInvMass);

      for (unsigned int i = 0; i < sizeInvMass; i++) {
        eigenvaluesInvariantMass(i) = eValIM[i];
        for (unsigned int j = 0; j < sizeMom; j++)
          eigenvectorsInvariantMass(i, j) = eVecIM[i * sizeMom + j];
        meanInvariantMass(i) = meanIM[i];
      }
    }

    /// Momentum parameters
    /// Generation-wise distributions
    Distributions momentumDistributions;
    /// Eigenvalues
    Acts::ActsDynamicVector eigenvaluesMomentum;
    /// Eigenvectors
    Acts::ActsDynamicMatrix eigenvectorsMomentum;
    /// Mean in eigenspace
    Acts::ActsDynamicVector meanMomentum;

    /// Invariant mass parameters
    /// Generation-wise distributions
    Distributions invariantMassDistributions;
    /// Eigenvalues
    Acts::ActsDynamicVector eigenvaluesInvariantMass;
    /// Eigenvectors
    Acts::ActsDynamicMatrix eigenvectorsInvariantMass;
    /// Mean in eigenspace
    Acts::ActsDynamicVector meanInvariantMass;
  };

  /// Initial momentum
  float momentum;
  /// Probability for soft nuclear interacion
  float softInteractionProbability;
  /// PDG ID based branching probabilities
  PdgMap pdgMap;
  /// Probability for nuclear interaction
  CumulativeDistribution nuclearInteractionProbability;
  /// Multiplicity in soft interactions
  CumulativeDistribution softMultiplicity;
  /// Multiplicity in hard interactions
  CumulativeDistribution hardMultiplicity;
  /// Kinematic distributions in soft interactions
  std::vector<ParametersWithFixedMultiplicity> softKinematicParameters;
  /// Kinematic distributions in hard interactions
  std::vector<ParametersWithFixedMultiplicity> hardKinematicParameters;
};

/// Parametrisation of a single particle
using Parametrisation = std::vector<std::pair<float, Parameters>>;
/// Parametrisation of multiple particles
using MultiParticleParametrisation =
    std::vector<std::pair<int, Parametrisation>>;
}  // namespace detail
}  // namespace ActsFatras
