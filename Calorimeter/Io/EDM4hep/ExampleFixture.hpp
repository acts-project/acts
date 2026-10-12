// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"

#include <cstdint>
#include <stdexcept>
#include <tuple>
#include <utility>

#include <edm4hep/CaloHitContributionCollection.h>
#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/SimCalorimeterHitCollection.h>
#include <podio/Frame.h>

namespace ActsExamples::CalorimeterFixture {

constexpr std::uint64_t highCellId = 0xfedcba9876543210ULL;

inline Acts::Vector3 cellCentre(std::uint64_t id) {
  if (id == highCellId) {
    return Acts::Vector3{1500, 10, 5} * Acts::UnitConstants::mm;
  }
  if (id == 42) {
    return Acts::Vector3{-1500, -10, -5} * Acts::UnitConstants::mm;
  }
  if (id == highCellId + 1) {
    return Acts::Vector3{1500, 20, 5} * Acts::UnitConstants::mm;
  }
  if (id == 7) {
    return Acts::Vector3{1500, 100, 5} * Acts::UnitConstants::mm;
  }
  throw std::invalid_argument("Unknown synthetic calorimeter cell");
}

/// Genuine EDM4hep collections in a synthetic fixture, without detector
/// transport.
inline podio::Frame makeFrame(bool extended = false) {
  edm4hep::MCParticleCollection particles;
  auto particle = particles.create();
  particle.setPDG(211);
  edm4hep::CaloHitContributionCollection contributions;
  edm4hep::SimCalorimeterHitCollection hits;
  auto first = hits.create();
  first.setCellID(highCellId);
  first.setEnergy(0.53f);
  first.setPosition({1, 2, 3});
  for (const auto [energy, time] :
       {std::pair{0.01f, 1.f}, std::pair{0.02f, 4.f}, std::pair{0.5f, 100.f}}) {
    auto contribution = contributions.create();
    contribution.setEnergy(energy);
    contribution.setTime(time);
    contribution.setParticle(particle);
    contribution.setStepPosition({4, 5, 6});
    first.addToContributions(contribution);
  }
  auto second = hits.create();
  second.setCellID(highCellId);
  second.setEnergy(0.01f);
  auto secondContribution = contributions.create();
  secondContribution.setEnergy(0.01f);
  secondContribution.setTime(7);
  secondContribution.setParticle(particle);
  second.addToContributions(secondContribution);
  auto third = hits.create();
  third.setCellID(42);
  third.setEnergy(0.02f);
  auto thirdContribution = contributions.create();
  thirdContribution.setEnergy(0.02f);
  thirdContribution.setTime(2);
  thirdContribution.setParticle(particle);
  third.addToContributions(thirdContribution);
  if (extended) {
    for (const auto [id, energy, time] :
         {std::tuple{highCellId + 1, 0.02f, 2.f},
          std::tuple{std::uint64_t{7}, 0.03f, 6.f}}) {
      auto hit = hits.create();
      hit.setCellID(id);
      hit.setEnergy(energy);
      auto contribution = contributions.create();
      contribution.setEnergy(energy);
      contribution.setTime(time);
      contribution.setParticle(particle);
      hit.addToContributions(contribution);
    }
  }
  podio::Frame frame;
  frame.put(std::move(particles), "MCParticles");
  frame.put(std::move(contributions), "CaloContributions");
  frame.put(std::move(hits), "SimCaloHits");
  return frame;
}

}  // namespace ActsExamples::CalorimeterFixture
