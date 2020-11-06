// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/HepMC/EventExtraction.hpp"

#include <unordered_map>
#include <TH1F.h>
#include <TVectorF.h>

namespace NuclearInteractionParametrisation {

struct EventFraction {
  EventFraction() = default;

  EventFraction(ActsExamples::SimParticle initPart,
                std::vector<ActsExamples::SimParticle> finalPart)
      : initialParticle(std::move(initPart)),
        finalParticles(std::move(finalPart)) {}

  ActsExamples::SimParticle initialParticle;
  std::vector<ActsExamples::SimParticle> finalParticles;

  bool soft = false;
  unsigned int multiplicity = 0;
};

static constexpr uint32_t s_MaxValue = UINT32_MAX;
using EventCollection = std::vector<EventFraction>;
using EventProperties = std::vector<std::vector<float>>;
using ProbabilityDistributions = std::vector<TH1F*>;
using CumulativeDistribution = TH1F*;
template <unsigned int length_t>
using Vector = Acts::ActsVectorF<length_t>;
template <unsigned int length_t>
using Matrix = Acts::ActsSymMatrixF<length_t>;
template <unsigned int length_t>
using EigenspaceComponents =
    std::tuple<Vector<length_t>, Matrix<length_t>, Vector<length_t>>;
template <unsigned int multiplicity_t>
using Parametrisation = std::pair<EigenspaceComponents<multiplicity_t + 1>,
                                  std::vector<CumulativeDistribution>>;

/// @brief This method scales the final state momenta by the initial momentum
///
/// @param [in] events The event storage
/// @param [in] multiplicity The target multiplicity
/// @param [in] soft Decision whether soft interactions should be considered
///
/// @return Vector containing the scaled final states
EventProperties prepateMomenta(const EventCollection& events,
                               unsigned int multiplicity, bool soft);

/// @brief This method builds the distributions of each generation of a final
/// state
///
/// @param [in] events The storage of a property with fixed multiplicity and
/// interaction type
/// @param [in] nBins Number of bins in the histogram
///
/// @return Vector containing the distributions for each generation
ProbabilityDistributions buildMomPerMult(const EventProperties& events,
                                         unsigned int nBins);

/// @brief Transform probability distribution entries to standard normal
/// distribution entries
///
/// @param [in] histos The probability distributions
/// @param [in] events The storage of a property with fixed multiplicity and
/// interaction type
///
/// @return The transformed values of @p events
EventProperties convertEventToGaussian(const ProbabilityDistributions& histos,
                                       const EventProperties& events);

/// @brief Calculate the mean and the covariance matrix of a final state
/// property
///
/// @tparam multiplicity_t The final state multiplicity + 1
/// @param [in] events The storage of a property with fixed multiplicity and
/// interaction type
///
/// @return Pair containing the mean and the covariance matrix
template <unsigned int multiplicity_t>
std::pair<Vector<multiplicity_t>, Matrix<multiplicity_t>>
calculateMeanAndCovariance(const EventProperties& events) {
  // Calculate the mean
  Vector<multiplicity_t> mean = Vector<multiplicity_t>::Zero();
  for (const std::vector<float>& event : events)
    for (unsigned int j = 0; j < multiplicity_t; j++)
      mean[j] += event[j];
  mean /= (float)events.size();

  // Calculate the covariance matrix
  Matrix<multiplicity_t> covariance = Matrix<multiplicity_t>::Zero();
  for (unsigned int i = 0; i < multiplicity_t; i++)
    for (unsigned int j = 0; j < multiplicity_t; j++)
      for (unsigned int k = 0; k < events.size(); k++)
        covariance(i, j) += (events[k][i] - mean[i]) * (events[k][j] - mean[j]);
  covariance /= (float)events.size();

  return std::make_pair(mean, covariance);
}

/// @brief Calculate the eigenvalues, eigenvectors and the mean in eigenspace
///
/// @tparam multiplicity_t The final state multiplicity + 1
/// @param [in] mean The mean of the normal distribution
/// @param [in] covariance The covariance matrix of the normal distribution
///
/// @return Tuple containing the eigenvalues, eigenvectors and the mean in
/// eigenspace
template <unsigned int multiplicity_t>
EigenspaceComponents<multiplicity_t> calculateEigenspace(
    const Vector<multiplicity_t>& mean,
    const Matrix<multiplicity_t>& covariance) {
  // Calculate eigenvalues and eigenvectors
  Eigen::EigenSolver<Matrix<multiplicity_t>> es(covariance);
  Vector<multiplicity_t> eigenvalues = es.eigenvalues().real();
  Matrix<multiplicity_t> eigenvectors = es.eigenvectors().real();
  // Transform the mean vector into eigenspace
  Vector<multiplicity_t> meanEigenspace = eigenvectors * mean;

  return std::make_tuple(eigenvalues, eigenvectors, meanEigenspace);
}

/// @brief This function calculates all components required for simulating final
/// state momenta
///
/// @tparam multiplicity_t The final state multiplicity that will be considered
/// @param [in] events The event storage
/// @param [in] soft Decision whether soft interactions should be considered
/// @param [in] nBins The number of bins in the histograms
///
/// @return Pair storing all components
template <unsigned int multiplicity_t>
Parametrisation<multiplicity_t> buildMomentumParameters(
    const EventCollection& events, bool soft, unsigned int nBins) {
  // Strip off data
  auto momenta = prepateMomenta(events, multiplicity_t, soft);

  // Build histos
  ProbabilityDistributions histos = buildMomPerMult(momenta, nBins);

  // Build normal distribution
  auto momentaGaussian = convertEventToGaussian(histos, momenta);
  auto meanAndCovariance =
      calculateMeanAndCovariance<multiplicity_t + 1>(momentaGaussian);
  // Calculate the transformation into the eigenspace of the covariance matrix
  EigenspaceComponents<multiplicity_t + 1> eigenspaceElements =
      calculateEigenspace<multiplicity_t + 1>(meanAndCovariance.first,
                                              meanAndCovariance.second);
  // Calculate the the cumulative distributions
  return std::make_pair(eigenspaceElements, histos);
}

/// @brief This method calculates the final state particles invariant masses
///
/// @param [in] events The event storage
/// @param [in] multiplicity The target multiplicity
/// @param [in] soft Decision whether soft interactions should be considered
///
/// @return Vector containing the final state invariant masses
EventProperties prepateInvariantMasses(const EventCollection& events,
                                       unsigned int multiplicity, bool soft);

/// @brief This function calculates all components required for simulating final
/// state invariant masses
///
/// @tparam multiplicity_t The final state multiplicity that will be considered
/// @param [in] events The event storage
/// @param [in] soft Decision whether soft interactions should be considered
/// @param [in] nBins The number of bins in the histograms
///
/// @return Pair storing all components
template <unsigned int multiplicity_t>
Parametrisation<multiplicity_t> buildInvariantMassParameters(
    const EventCollection& events, bool soft, unsigned int nBins) {
  // Strip off data
  auto invariantMasses = prepateInvariantMasses(events, multiplicity_t, soft);

  // Build histos
  ProbabilityDistributions histos = buildMomPerMult(invariantMasses, nBins);

  // Build normal distribution
  auto invariantMassesGaussian =
      convertEventToGaussian(histos, invariantMasses);
  auto meanAndCovariance =
      calculateMeanAndCovariance<multiplicity_t + 1>(invariantMassesGaussian);
  // Calculate the transformation into the eigenspace of the covariance matrix
  EigenspaceComponents<multiplicity_t + 1> eigenspaceElements =
      calculateEigenspace<multiplicity_t + 1>(meanAndCovariance.first,
                                              meanAndCovariance.second);
  // Calculate the the cumulative distributions
  return std::make_pair(eigenspaceElements, histos);
}

/// @brief This method evaluates the cumulative probabilities for a given
/// particle type to produce a particle type
///
/// @param [in] events The event storage
///
/// @return Map containing the cumulative branching probabilities
std::unordered_map<int, std::unordered_map<int, float>>
cumulativePDGprobability(const EventCollection& events);

/// @brief Evaluates the cumulative probabilities for a certain multiplicity of
/// in a soft or hard process
///
/// @param [in] events The event storage
///
/// @return Pair containing the distribution for soft and hard processes
std::pair<CumulativeDistribution, CumulativeDistribution>
cumulativeMultiplicityProbability(const EventCollection& events);

/// @brief This method evaluates the probability that a nuclear interaction is a
/// soft interaction
///
/// @param [in] events The event storage
///
/// @return The probability for soft interactions
TVectorF softProbability(const EventCollection& events);

/// @brief This method calculates the cumulative probability for a nuclear
/// interaction as a function of L0
///
/// @param [in] events The event storage
/// @param [in] interactionProbabilityBins Number of bins used for the histogram
/// @note The number of events is used for the normalisation of the
/// distribution. Hence the result is not necessarily normalised to 1. This
/// allows to sample whether a nuclear interaction occurs up to a certain value
/// of L0 or not at all.
///
/// @return The cumulative distribution for the nuclear interaction
CumulativeDistribution cumulativeNuclearInteractionProbability(
    const EventCollection& events, unsigned int interactionProbabilityBins);
}  // namespace NuclearInteractionParametrisation