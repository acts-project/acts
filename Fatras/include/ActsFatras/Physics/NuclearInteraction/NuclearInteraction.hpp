// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Common.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"
#include "ActsFatras/Physics/NuclearInteraction/Parameters.hpp"

#include <any>
#include <cmath>
#include <iterator>
#include <limits>
#include <random>
#include <vector>

namespace ActsFatras {

/// @brief This class provides a parametrised nuclear interaction. The thereby
/// required parametrisation needs to be set and is not provided by default.
///
/// @note This class differs between two different processes labelled as nuclear
/// interaction. Either the initial particle survives (soft) or it gets
/// destroyed (hard) by this process.
struct NuclearInteraction {
  /// The storage of the parameterisation
  detail::MultiParticleParametrisation multiParticleParameterisation;
  /// The number of trials to match momenta and inveriant masses
  unsigned int nMatchingTrials = std::numeric_limits<unsigned int>::max();

  /// @brief Main call operator
  ///
  /// @tparam generator_t The random number generator type
  /// @param [in, out] generator The random number generator
  /// @param [in] slab The material
  /// @param [in, out] particle The ingoing particle
  ///
  /// @return Vector containing the produced secondaries
  template <typename generator_t>
  std::vector<Particle> operator()(generator_t& generator,
                                   const Acts::MaterialSlab& /*slab*/,
                                   Particle& particle) const {
    // Fast exit: No paramtrization provided
    if (multiParticleParameterisation.empty()) {
      return {};
    }

    // Find the parametrisation that corresponds to the particle type
    for (const auto& particleParametrisation : multiParticleParameterisation) {
      if (particleParametrisation.first == particle.pdg()) {
        std::uniform_real_distribution<double> uniformDistribution{0., 1.};

        // Get the parameters
        const detail::Parametrisation& singleParticleParametrisation =
            particleParametrisation.second;
        const detail::Parameters& parametrisation = findParameters(
            uniformDistribution(generator), singleParticleParametrisation,
            particle.absoluteMomentum());

        // Set the L0 limit if not done already
        if (particle.pathLimitL0() ==
            std::numeric_limits<Particle::Scalar>::max()) {
          const auto& distribution =
              parametrisation.nuclearInteractionProbability;
          particle.setMaterialLimits(
              particle.pathLimitX0(),
              sampleContinuousValues(uniformDistribution(generator),
                                     distribution));
        }

        // Test whether enough material was passed for a nuclear interaction
        if (particle.pathInL0() >= particle.pathLimitL0()) {
          std::normal_distribution<double> normalDistribution{0., 1.};
          // Dice the interaction type
          const bool interactSoft =
              softInteraction(normalDistribution(generator),
                              parametrisation.softInteractionProbability);

          // Get the final state multiplicity
          const unsigned int multiplicity = finalStateMultiplicity(
              uniformDistribution(generator),
              interactSoft ? parametrisation.softMultiplicity
                              : parametrisation.hardMultiplicity);

          // Get the particle types
          const std::vector<int> pdgIds =
              samplePdgIds(generator, parametrisation.pdgMap, multiplicity,
                           particle.pdg(), interactSoft);

          // Get the multiplicity
          const std::vector<
              detail::Parameters::ParametersWithFixedMultiplicity>&
              parametrisationOfType =
                  interactSoft ? parametrisation.softKinematicParameters
                                  : parametrisation.hardKinematicParameters;
          const detail::Parameters::ParametersWithFixedMultiplicity&
              parametrisationOfMultiplicity =
                  parametrisationOfType[multiplicity];

          // Get the kinematics
          const auto kinematics = sampleKinematics(
              generator, parametrisationOfMultiplicity,
              parametrisation.momentum);

          // Construct the particles
          const auto particles = convertParametersToParticles(
              generator, pdgIds, kinematics.first, kinematics.second, particle,
              interactSoft);

          // Kill the particle in a hard process
          if (!interactSoft)
            particle.setAbsoluteMomentum(0);

          return particles;
        }
      } else {
        // Fast exit if no parametrisation for the particle was provided
        return {};
      }
    }
    return {};
  }

 private:
  /// @brief Retrieves the parametrisation for the particle
  ///
  /// @param [in] rnd A random number
  /// @param [in] parametrisation The storage of parametrisations
  /// @param [in] particleMomentum The particles momentum
  ///
  /// @return The parametrisation
  const detail::Parameters& findParameters(
      double rnd, const detail::Parametrisation& parametrisation,
      float particleMomentum) const;

  /// @brief Estimates the interaction type
  ///
  /// @param [in] rnd Random number
  /// @param [in] probability The probability for a soft interaction
  ///
  /// @return True if a soft interaction occurs
  inline bool softInteraction(double rnd, float probability) const {
    return rnd <= probability;
  }

  /// @brief Evaluates the multiplicity of the final state
  ///
  /// @param [in] rnd Random number
  /// @param [in] distribution The multiplicity distribution
  ///
  /// @return The final state multiplicity
  unsigned int finalStateMultiplicity(
      double rnd,
      const detail::Parameters::CumulativeDistribution& distribution) const;

  /// @brief Evaluates the final state PDG IDs
  ///
  /// @tparam generator_t The random number generator type
  /// @param [in, out] generator The random number generator
  /// @param [in] pdgMap The branching probability map
  /// @param [in] multiplicity The final state multiplicity
  /// @param [in] particlePdg The PDG ID of the initial particle
  /// @param [in] soft Treat it as soft or hard nuclear interaction
  ///
  /// @return Vector containing the PDG IDs
  template <typename generator_t>
  std::vector<int> samplePdgIds(generator_t& generator,
                                const detail::Parameters::PdgMap& pdgMap,
                                unsigned int multiplicity, int particlePdg,
                                bool soft) const;

  /// @brief Evaluates the final state invariant masses
  ///
  /// @tparam generator_t The random number generator type
  /// @param [in, out] generator The random number generator
  /// @param [in] kinematicParameters Parametrisation of kinematic properties
  ///
  /// @return Vector containing the invariant masses
  template <typename generator_t>
  Acts::ActsDynamicVector sampleInvariantMasses(
      generator_t& generator,
      const detail::Parameters::ParametersWithFixedMultiplicity&
          parametrisation) const;

  /// @brief Evaluates the final state momenta
  ///
  /// @tparam generator_t The random number generator type
  /// @param [in, out] The random number generator
  /// @param [in] kinematicParameters Parametrisation of kinematic properties
  ///
  /// @return Vector containing the momenta
  template <typename generator_t>
  Acts::ActsDynamicVector sampleMomenta(
      generator_t& generator,
      const detail::Parameters::ParametersWithFixedMultiplicity&
          parametrisation,
      float initialMomentum) const;

  /// @brief Tests whether the final state momenta and invariant masses are
  /// matching to each other to allow the evaluation of particle directions.
  ///
  /// @param [in] momenta The final state momenta
  /// @param [in] invariantMasses The final state invariant masses
  /// @param [in] initialMomentum The momentum of the initial particle
  ///
  /// @return Decision whether the parameters can be matched to each other or
  /// not.
  bool match(const Acts::ActsDynamicVector& momenta,
             const Acts::ActsDynamicVector& invariantMasses,
             float initialMomentum) const;

  /// @brief This method samples the kinematics of the final state particles
  ///
  /// @tparam generator_t The random number generator type
  /// @param [in, out] The random number generator
  /// @param [in] parameters The parametrisation
  /// @param [in] momentum The momentum of the parametrisation
  ///
  /// @return The final state momenta and invariant masses
  template <typename generator_t>
  std::pair<Acts::ActsDynamicVector, Acts::ActsDynamicVector> sampleKinematics(
      generator_t& generator,
      const detail::Parameters::ParametersWithFixedMultiplicity& parameters,
      float momentum) const;

  /// @brief Converts relative angles to absolute angles wrt the global
  /// coordinate system.
  /// @note It is assumed that the angles of the first particle are provided in
  /// the context of the global coordinate system whereas the angles of the
  /// second particle are provided relatively to the first particle.
  ///
  /// @param [in] phi1 The azimuthal angle of the first particle
  /// @param [in] theta1 The polar angle of the first particle
  /// @param [in] phi2 The azimuthal angle of the second particle
  /// @param [in] theta2 The polar angle of the second particle
  ///
  /// @return Azimuthal and polar angle of the second particle in the global
  /// coordinate system
  std::pair<ActsFatras::Particle::Scalar, ActsFatras::Particle::Scalar>
  globalAngle(ActsFatras::Particle::Scalar phi1,
              ActsFatras::Particle::Scalar theta1, float phi2,
              float theta2) const;

  /// @brief Converter from sampled numbers to a vector of particles
  ///
  /// @tparam generator_t The random number generator type
  /// @param [in, out] The random number generator
  /// @param [in] pdgId The PDG IDs
  /// @param [in] momentum The momenta
  /// @param [in] invariantMass The invariant masses
  /// @param [in] initialParticle The initial particle
  /// @param [in] soft Treat it as soft or hard nuclear interaction
  ///
  /// @return Vector containing the final state particles
  template <typename generator_t>
  std::vector<Particle> convertParametersToParticles(
      generator_t& generator, const std::vector<int>& pdgId,
      const Acts::ActsDynamicVector& momenta,
      const Acts::ActsDynamicVector& invariantMasses, Particle& initialParticle,
      bool soft) const;

  /// @brief This function performs an inverse sampling to provide a discrete
  /// value from a distribution.
  ///
  /// @param [in] rnd A random number in [0,1]
  /// @param [in] distribution The distribution to sample from
  ///
  /// @return The sampled value
  unsigned int sampleDiscreteValues(
      double rnd,
      const detail::Parameters::CumulativeDistribution& distribution) const;

  /// @brief This function performs an inverse sampling to provide a continuous
  /// value from a distribition.
  ///
  /// @param [in] rnd A random number in [0,1]
  /// @param [in] distribution The distribution to sample from
  /// @param [in] interpolate Flag to steer whether an interpolation between
  /// neighbouring bins should be performed instead of a bin lookup
  ///
  /// @return The sampled value
  double sampleContinuousValues(
      double rnd,
      const detail::Parameters::CumulativeDistribution& distribution,
      bool interpolate = false) const;
};

template <typename generator_t>
std::vector<int> NuclearInteraction::samplePdgIds(
    generator_t& generator, const detail::Parameters::PdgMap& pdgMap,
    unsigned int multiplicity, int particlePdg, bool soft) const {
  // Fast exit in case of no final state particles
  if (multiplicity == 0)
    return {};

  // The final state PDG IDs
  std::vector<int> pdgIds;
  pdgIds.reserve(multiplicity);

  std::uniform_real_distribution<float> uniformDistribution{0., 1.};

  // Find the producers probability distribution
  auto citProducer = pdgMap.cbegin();
  while (citProducer->first != particlePdg && citProducer != pdgMap.end())
    citProducer++;

  const std::vector<std::pair<int, float>>& mapInitial = citProducer->second;

  // Set the first particle depending on the interaction type
  if (soft)
    // Store the initial particle if the interaction is soft
    pdgIds.push_back(particlePdg);
  else {
    // Otherwise dice the particle
    const float rndInitial = uniformDistribution(generator);
    pdgIds.push_back(
        std::lower_bound(mapInitial.begin(), mapInitial.end(), rndInitial,
                         [](const std::pair<int, float>& element,
                            float random) { return element.second < random; })
            ->second);
  }

  // Set the remaining particles
  for (unsigned int i = 1; i < multiplicity; i++) {
    // Find the producers probability distribution from the last produced
    // particle
    citProducer = pdgMap.cbegin();
    while (citProducer->first != pdgIds[i - 1] && citProducer != pdgMap.end())
      citProducer++;

    // Set the next particle
    const std::vector<std::pair<int, float>>& map = citProducer->second;
    const float rnd = uniformDistribution(generator);
    pdgIds.push_back(
        std::lower_bound(map.begin(), map.end(), rnd,
                         [](const std::pair<int, float>& element,
                            float random) { return element.second < random; })
            ->second);
  }
  return pdgIds;
}

template <typename generator_t>
Acts::ActsDynamicVector NuclearInteraction::sampleInvariantMasses(
    generator_t& generator,
    const detail::Parameters::ParametersWithFixedMultiplicity& parametrisation)
    const {
  // The resulting vector
  Acts::ActsDynamicVector parameters;
  const unsigned int size = parametrisation.eigenvaluesInvariantMass.size();
  parameters.resize(size);

  // Sample in the eigenspace
  for (unsigned int i = 0; i < size; i++) {
    float variance = parametrisation.eigenvaluesInvariantMass[i];
    std::normal_distribution<Acts::ActsScalar> dist{parametrisation.meanInvariantMass[i],
                                         sqrtf(variance)};
    parameters[i] = dist(generator);
  }
  // Transform to multivariate normal distribution
  parameters = parametrisation.eigenvectorsInvariantMass * parameters;

  // Perform the inverse sampling from the distributions
  for (unsigned int i = 0; i < size; i++) {
    const double cdf = (std::erff(parameters[i]) + 1) * 0.5;
    parameters[i] = sampleContinuousValues(
        cdf, parametrisation.invariantMassDistributions[i]);
  }
  return parameters;
}

template <typename generator_t>
Acts::ActsDynamicVector NuclearInteraction::sampleMomenta(
    generator_t& generator,
    const detail::Parameters::ParametersWithFixedMultiplicity& parametrisation,
    float initialMomentum) const {
  // The resulting vector
  Acts::ActsDynamicVector parameters;
  const unsigned int size = parametrisation.eigenvaluesMomentum.size();
  parameters.resize(size);

  // Sample in the eigenspace
  for (unsigned int i = 0; i < size + 1; i++) {
    float variance = parametrisation.eigenvaluesMomentum[i];
    std::normal_distribution<Acts::ActsScalar> dist{parametrisation.meanMomentum[i],
                                         sqrtf(variance)};
    parameters[i] = dist(generator);
  }
  // Transform to multivariate normal distribution
  parameters = parametrisation.eigenvectorsMomentum * parameters;

  // Perform the inverse sampling from the distributions
  for (unsigned int i = 0; i < size + 1; i++) {
    const float cdf = (std::erff(parameters[i]) + 1) * 0.5;
    parameters[i] =
        sampleContinuousValues(cdf, parametrisation.momentumDistributions[i]);
  }

  // Scale the momenta
  Acts::ActsDynamicVector momenta = parameters.template head(size);
  const float sum = momenta.sum();
  const float scale = parameters.template tail<1>()(0, 0) / sum;
  momenta *= scale * initialMomentum;

  return momenta;
}

template <typename generator_t>
std::pair<Acts::ActsDynamicVector, Acts::ActsDynamicVector>
NuclearInteraction::sampleKinematics(
    generator_t& generator,
    const detail::Parameters::ParametersWithFixedMultiplicity& parameters,
    float momentum) const {
  unsigned int trials = 0;
  auto invariantMasses = sampleInvariantMasses(generator, parameters);
  auto momenta = sampleMomenta(generator, parameters, momentum);
  // Repeat momentum evaluation until the parameters match
  while (!match(momenta, invariantMasses, momentum)) {
    momenta = sampleMomenta(generator, parameters, momentum);
    // Re-sampole invariant masses if no fitting momenta were found
    if (trials % nMatchingTrials == 0)
      invariantMasses = sampleInvariantMasses(generator, parameters);
  }
  return std::make_pair(momenta, invariantMasses);
}

template <typename generator_t>
std::vector<Particle> NuclearInteraction::convertParametersToParticles(
    generator_t& generator, const std::vector<int>& pdgId,
    const Acts::ActsDynamicVector& momenta,
    const Acts::ActsDynamicVector& invariantMasses, Particle& initialParticle,
    bool soft) const {
  std::uniform_real_distribution<double> uniformDistribution{0., 1.};
  const auto initialMomentum = initialParticle.absoluteMomentum();
  const auto& initialDirection = initialParticle.unitDirection();
  const double phi = Acts::VectorHelpers::phi(initialDirection);
  const double theta = Acts::VectorHelpers::theta(initialDirection);
  const unsigned int size = momenta.size();

  std::vector<Particle> result;
  result.reserve(size);

  // Build the particles
  for (unsigned int i = 0; i < size; i++) {
    const float momentum = momenta[i];
    const float invariantMass = invariantMasses[i];
    const float p1p2 = 2. * momentum * initialMomentum;
    const float costheta = 1. - invariantMass * invariantMass / p1p2;

    const auto phiTheta =
        globalAngle(phi, theta, uniformDistribution(generator) * 2. * M_PI,
                    std::acos(costheta));
    const auto direction =
        Acts::makeDirectionUnitFromPhiTheta(phiTheta.first, phiTheta.second);
    Particle p = Particle(Barcode(), static_cast<Acts::PdgParticle>(pdgId[i]));
    p.setProcess(ProcessType::eNuclearInteraction)
        .setPosition4(initialParticle.fourPosition())
        .setAbsoluteMomentum(momentum)
        .setDirection(direction);

    // Search for the parametrisation of the particle
    auto cit = multiParticleParameterisation.begin();
    while (cit->first != p.pdg() && cit != multiParticleParameterisation.end())
      cit++;

    if (cit != multiParticleParameterisation.end()) {
      // Assign a path limit in L0 to the particle
      const auto& distribution =
          findParameters(uniformDistribution(generator), cit->second,
                         p.absoluteMomentum())
              .nuclearInteractionProbability;
      p.setMaterialLimits(
          p.pathLimitX0(),
          sampleContinuousValues(uniformDistribution(generator), distribution));
    }

    // Store the particle
    if (i == 0 && soft)
      initialParticle = p;
    else
      result.push_back(std::move(p));
  }

  return result;
}
}  // namespace ActsFatras
