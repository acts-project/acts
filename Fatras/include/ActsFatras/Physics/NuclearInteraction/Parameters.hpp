// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"
#include <map>
#include <unordered_map>
#include <vector>

namespace ActsFatras {
namespace detail {

/// @brief Data storage of the parametrised nuclear interaction
struct Parameters {
  using CumulativeDistribution =
      std::pair<std::vector<float>, std::vector<uint32_t>>;
  using Distributions = std::vector<CumulativeDistribution>;
  using PdgMap = std::unordered_map<int, std::unordered_map<int, float>>;

  /// @brief Nested struct for the storage of the kinematic parametrisations for
  /// a given final state multiplicity
  struct ParametersWithFixedMultiplicity {
    ParametersWithFixedMultiplicity(
        Distributions&& momenta, std::vector<float>& eValMom,
        std::vector<float>& eVecMom, std::vector<float>& meanMom,
        Distributions&& invariantMasses, std::vector<float>& eValIM,
        std::vector<float>& eVecIM, std::vector<float>& meanIM)
        : momentumDistributions(momenta),
          eigenvaluesMomentum(eValMom.data()),
          eigenvectorsMomentum(eVecMom.data()),
          meanMomentum(meanMom.data()),
          invariantMassDistributions(invariantMasses),
          eigenvaluesInvariantMass(eValIM.data()),
          eigenvectorsInvariantMass(eVecIM.data()),
          meanInvariantMass(meanIM.data()) {}

    /// Momentum parameters
    /// Generation-wise distributions
    Distributions momentumDistributions;
    /// Eigenvalues
    Acts::ActsVectorXf eigenvaluesMomentum;
    /// Eigenvectors
    Acts::ActsMatrixXf eigenvectorsMomentum;
    /// Mean in eigenspace
    Acts::ActsVectorXf meanMomentum;

    /// Invariant mass parameters
    /// Generation-wise distributions
    Distributions invariantMassDistributions;
    /// Eigenvalues
    Acts::ActsVectorXf eigenvaluesInvariantMass;
    /// Eigenvectors
    Acts::ActsMatrixXf eigenvectorsInvariantMass;
    /// Mean in eigenspace
    Acts::ActsVectorXf meanInvariantMass;
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
using Parametrisation = std::map<float, Parameters>;
/// Parametrisation of multiple particles
using MultiParticleParametrisation = std::unordered_map<int, Parametrisation>;
}  // namespace detail
}  // namespace ActsFatras
