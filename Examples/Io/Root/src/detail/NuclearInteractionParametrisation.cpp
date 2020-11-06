// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/detail/NuclearInteractionParametrisation.hpp"

#include <limits>
#include <Eigen/Eigenvalues>
#include <TMath.h>

namespace NuclearInteractionParametrisation {
namespace {

/// @brief Evaluate the location in a standard normal distribution for a value
/// from a probability distribution
///
/// @param [in] histo The probability distribution
/// @param [in] mom The abscissa value in @p histo
///
/// @return The location in a standard normal distribution
float gaussianValue(TH1F const* histo, const float mom) {
  // Get the cumulative probability distribution
  TH1F* normalised = (TH1F*)histo->DrawNormalized();
  TH1F* cumulative = (TH1F*)normalised->GetCumulative();
  // Find the cumulative probability
  const float binContent = cumulative->GetBinContent(cumulative->FindBin(mom));
  // Transform the probability to an entry in a standard normal distribution
  const float value = TMath::ErfInverse(2. * binContent - 1.);

  delete (normalised);
  delete (cumulative);
  return value;
}
}  // namespace

namespace {
float invariantMass(const ActsExamples::SimParticle::Vector4& fourVector1,
                    const ActsExamples::SimParticle::Vector4& fourVector2) {
  ActsExamples::SimParticle::Vector4 sum = fourVector1 + fourVector2;
  const ActsExamples::SimParticle::Scalar energy = sum[Acts::eEnergy];
  ActsExamples::SimParticle::Vector3 momentum =
      sum.template segment<3>(Acts::eMom0);
  return std::sqrt(energy * energy - momentum.norm());
}
}  // namespace

EventProperties prepateMomenta(const EventCollection& events,
                               unsigned int multiplicity,
                               bool soft)  // TODO: build enum instead of bool
{
  EventProperties result;
  // Loop over all events
  for (const EventFraction& event : events) {
    // Test the multiplicity and type of the event
    if (event.multiplicity == multiplicity && event.soft == soft) {
      const float initialMomentum = event.initialParticle.absMomentum();
      float sum = 0.;
      std::vector<float> momenta;
      momenta.reserve(multiplicity);
      // Fill the vector with the scaled momenta
      for (const ActsExamples::SimParticle& p : event.finalParticles) {
        sum += p.absMomentum();
        momenta.push_back(p.absMomentum() / initialMomentum);
      }
      // Add the scaled sum of momenta
      momenta.push_back(sum / initialMomentum);
      result.push_back(std::move(momenta));
    }
  }
  return result;
}

ProbabilityDistributions buildMomPerMult(const EventProperties& events,
                                         unsigned int nBins) {
  // Fast exit
  if (events.empty())
    return {};
  const unsigned int multMax = events[0].size();

  // Find the range of each histogram
  std::vector<float> min(multMax, std::numeric_limits<float>::max());
  std::vector<float> max(multMax, 0);
  for (const std::vector<float>& event : events)
    for (unsigned int i = 0; i < multMax; i++) {
      min[i] = std::min(event[i], min[i]);
      max[i] = std::max(event[i], max[i]);
    }

  // Evaluate the range of the histograms
  // This is used to avoid entries in over-/underflow bins
  std::vector<float> diff(multMax);
  for (unsigned int i = 0; i < multMax; i++)
    diff[i] = (max[i] - min[i]) * 0.1;

  // Build the histograms
  ProbabilityDistributions histos(multMax);
  for (unsigned int i = 0; i < multMax; i++) {
    histos[i] = new TH1F("", "", nBins, min[i] - diff[i], max[i] + diff[i]);
  }

  // Fill the histograms
  for (const std::vector<float>& event : events) {
    for (unsigned int i = 0; i < multMax; i++) {
      histos[i]->Fill(event[i]);
    }
  }
  return histos;
}

EventProperties convertEventToGaussian(const ProbabilityDistributions& histos,
                                       const EventProperties& events) {
  // Fast exit
  if (events.empty())
    return {};
  const unsigned int multMax = events[0].size();

  // Loop over the events
  EventProperties gaussianEvents;
  for (const std::vector<float>& event : events) {
    // Transform the properties in the events
    std::vector<float> gaussianEvent;
    for (unsigned int i = 0; i < multMax; i++) {
      gaussianEvent.push_back(gaussianValue(histos[i], event[i]));
    }
    // Store the transformed event
    gaussianEvents.push_back(gaussianEvent);
  }
  return gaussianEvents;
}

EventProperties prepateInvariantMasses(const EventCollection& events,
                                       unsigned int multiplicity, bool soft) {
  EventProperties result;
  // Loop over all events
  for (const EventFraction& event : events) {
    // Test the multiplicity and type of the event
    if (event.multiplicity == multiplicity && event.soft == soft) {
      const auto initialFourVector = event.initialParticle.momentum4();
      std::vector<float> invariantMasses;
      invariantMasses.reserve(multiplicity);
      // Fill the vector with the invariant masses
      for (const ActsExamples::SimParticle& p : event.finalParticles) {
        const auto fourVector = p.momentum4();
        invariantMasses.push_back(invariantMass(initialFourVector, fourVector));
      }
      result.push_back(invariantMasses);
    }
  }
  return result;
}

std::unordered_map<int, std::unordered_map<int, float>>
cumulativePDGprobability(const EventCollection& events) {
  std::unordered_map<int, std::unordered_map<int, float>> counter;
  std::unordered_map<int, float> totalSum;

  for (const EventFraction& event : events) {
    if (!event.soft) {
      counter[event.initialParticle.pdg()][event.finalParticles[0].pdg()]++;
      totalSum[event.initialParticle.pdg()]++;
    }
    for (unsigned int i = 1; i < event.multiplicity; i++) {
      counter[event.finalParticles[i - 1].pdg()]
             [event.finalParticles[i].pdg()]++;
      totalSum[event.finalParticles[i - 1].pdg()]++;
    }
  }

  // Build a cumulative distribution
  for (const auto& element : counter) {
    for (auto it1 = counter[element.first].begin();
         it1 != counter[element.first].end(); it1++) {
      // Add content to next bins
      for (auto it2 = std::next(it1, 1); it2 != counter[element.first].end();
           it2++) {
        it2->second += it1->second;
      }
      // Normalise the entry
      it1->second /= totalSum[element.first];
    }
  }
  return counter;
}

std::pair<CumulativeDistribution, CumulativeDistribution>
cumulativeMultiplicityProbability(const EventCollection& events) {
  // Find the range of both histogram
  unsigned int minSoft = std::numeric_limits<unsigned int>::max();
  unsigned int maxSoft = 0;
  unsigned int minHard = std::numeric_limits<unsigned int>::max();
  unsigned int maxHard = 0;
  for (const EventFraction& event : events) {
    if (event.soft) {
      minSoft = std::min(event.multiplicity, minSoft);
      maxSoft = std::max(event.multiplicity, maxSoft);
    } else {
      minHard = std::min(event.multiplicity, minHard);
      maxHard = std::max(event.multiplicity, maxHard);
    }
  }

  // Build and fill the histograms
  TH1F* softHisto = new TH1F("", "", maxSoft + 1, minSoft, maxSoft);
  TH1F* hardHisto = new TH1F("", "", maxHard + 1, minHard, maxHard);
  for (const EventFraction& event : events) {
    if (event.soft)
      softHisto->Fill(event.multiplicity);
    else
      hardHisto->Fill(event.multiplicity);
  }

  return std::make_pair(softHisto, hardHisto);
}

TVectorF softProbability(const EventCollection& events) {
  float counter = 0.;
  // Count the soft events
  for (const EventFraction& event : events)
    if (event.soft)
      counter++;

  TVectorF result(1);
  result[0] = counter / (float)events.size();
  return result;
}

CumulativeDistribution cumulativeNuclearInteractionProbability(
    const EventCollection& events, unsigned int interactionProbabilityBins) {
  // Find the limits of the histogram
  float min = std::numeric_limits<float>::max();
  float max = 0.;
  for (const EventFraction& event : events) {
    min = std::min((float)event.initialParticle.pathInL0(), min);
    max = std::max((float)event.initialParticle.pathInL0(), max);
  }

  // Fill the histogram
  TH1F* histo = new TH1F("", "", interactionProbabilityBins, min, max);
  for (const EventFraction& event : events)
    histo->Fill(event.initialParticle.pathInL0());

  // Build the distributions
  return histo;  // TODO: in this case the normalisation is not taking into
                 // account
}
}  // namespace NuclearInteractionParametrisation