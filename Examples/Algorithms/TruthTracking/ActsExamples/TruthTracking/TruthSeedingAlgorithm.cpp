// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/TruthSeedingAlgorithm.hpp"

#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Utilities/Range.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace ActsExamples {

TruthSeedingAlgorithm::TruthSeedingAlgorithm(Config cfg,
                                             Acts::Logging::Level lvl)
    : IAlgorithm("TruthSeedingAlgorithm", lvl), m_cfg(std::move(cfg)) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input truth particles collection");
  }
  if (m_cfg.inputParticleMeasurementsMap.empty()) {
    throw std::invalid_argument(
        "Missing input particle-measurements map collection");
  }
  if (m_cfg.inputSpacePoints.empty()) {
    throw std::invalid_argument("Missing seeds or space point collection");
  }

  for (const auto& spName : m_cfg.inputSpacePoints) {
    if (spName.empty()) {
      throw std::invalid_argument("Invalid space point input collection");
    }

    auto& handle = m_inputSpacePoints.emplace_back(
        std::make_unique<ReadDataHandle<SimSpacePointContainer>>(
            this,
            "InputSpacePoints#" + std::to_string(m_inputSpacePoints.size())));
    handle->initialize(spName);
  }

  if (m_cfg.outputParticles.empty()) {
    throw std::invalid_argument("Missing output particles collection");
  }
  if (m_cfg.outputSeeds.empty()) {
    throw std::invalid_argument("Missing seeds output collections");
  }
  if (m_cfg.outputProtoTracks.empty()) {
    throw std::invalid_argument("Missing proto tracks output collections");
  }

  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing input simulated hits collection");
  }

  if (m_cfg.inputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument(
        "Missing input simulated hits measurements map");
  }

  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputParticleMeasurementsMap.initialize(m_cfg.inputParticleMeasurementsMap);
  m_inputSimHits.initialize(m_cfg.inputSimHits);
  m_inputMeasurementSimHitsMap.initialize(m_cfg.inputMeasurementSimHitsMap);
  m_outputParticles.initialize(m_cfg.outputParticles);
  m_outputProtoTracks.initialize(m_cfg.outputProtoTracks);
  m_outputSeeds.initialize(m_cfg.outputSeeds);
  m_outputParticleHypotheses.maybeInitialize(m_cfg.outputParticleHypotheses);
}

ProcessCode TruthSeedingAlgorithm::execute(const AlgorithmContext& ctx) const {
  // prepare input collections
  const auto& particles = m_inputParticles(ctx);
  const auto& particleMeasurementsMap = m_inputParticleMeasurementsMap(ctx);
  const auto& simHits = m_inputSimHits(ctx);
  const auto& measurementSimHitsMap = m_inputMeasurementSimHitsMap(ctx);

  // construct the combined input container of space point pointers from all
  // configured input sources.
  // pre-compute the total size required so we only need to allocate once
  std::size_t nSpacePoints = 0;
  for (const auto& isp : m_inputSpacePoints) {
    nSpacePoints += (*isp)(ctx).size();
  }

  std::vector<const SimSpacePoint*> spacePointPtrs;
  spacePointPtrs.reserve(nSpacePoints);
  for (const auto& isp : m_inputSpacePoints) {
    for (const auto& spacePoint : (*isp)(ctx)) {
      // since the event store owns the space points, their pointers should be
      // stable and we do not need to create local copies.
      spacePointPtrs.push_back(&spacePoint);
    }
  }

  SimParticleContainer seededParticles;
  SimSeedContainer seeds;
  ProtoTrackContainer tracks;
  std::vector<Acts::ParticleHypothesis> particleHypotheses;

  seededParticles.reserve(particles.size());
  seeds.reserve(particles.size());
  tracks.reserve(particles.size());
  particleHypotheses.reserve(particles.size());

  std::unordered_map<Index, const SimSpacePoint*> spMap;

  for (const auto& spp : spacePointPtrs) {
    if (spp->sourceLinks().empty()) {
      ACTS_WARNING("Missing source link in space point");
      continue;
    }
    for (const auto& slink : spp->sourceLinks()) {
      const IndexSourceLink& islink = slink.get<IndexSourceLink>();
      spMap.emplace(islink.index(), spp);
    }
  }

  for (const auto& particle : particles) {
    // find the corresponding measurements for this particle
    const auto& measurements =
        makeRange(particleMeasurementsMap.equal_range(particle.particleId()));
    // fill measurement indices to create the prototrack
    ProtoTrack track;
    track.reserve(measurements.size());

    std::vector<std::pair<const SimHit*, Index>> hits;
    hits.reserve(measurements.size());

    for (const auto& [barcode, index] : measurements) {
      const auto simHitMapIt = measurementSimHitsMap.find(index);
      if (simHitMapIt == measurementSimHitsMap.end()) {
        ACTS_WARNING("No sim hit found for measurement index " << index);
        continue;
      }

      const auto simHitIt = simHits.nth(simHitMapIt->second);
      if (simHitIt == simHits.end()) {
        ACTS_WARNING("No sim hit found for index " << simHitMapIt->second);
        continue;
      }

      const auto& simHit = *simHitIt;

      hits.emplace_back(&simHit, index);
    }

    std::sort(hits.begin(), hits.end(), [](const auto& a, const auto& b) {
      return a.first->time() < b.first->time();
    });

    for (const auto& [hit, index] : hits) {
      track.push_back(index);
    }

    // The list of measurements and the initial start parameters
    if (track.size() < 3) {
      ACTS_WARNING("Particle " << particle << " has less than 3 measurements");
      continue;
    }
    // Space points on the prototrack
    std::vector<const SimSpacePoint*> spacePointsOnTrack;
    spacePointsOnTrack.reserve(track.size());
    // Loop over the measurement index on the prototrack to find the space
    // points
    for (const auto& measurementIndex : track) {
      auto it = spMap.find(measurementIndex);
      if (it != spMap.end()) {
        spacePointsOnTrack.push_back(it->second);
      }
    }
    // At least three space points are required
    if (spacePointsOnTrack.size() < 3) {
      continue;
    }

    // Loop over the found space points to find the seed with the maximum score.
    // The score is defined as the product of the deltaR of the bottom-middle
    // and middle-top space point pairs to favor separation between both pairs.
    bool seedFound = false;
    std::array<std::size_t, 3> bestSPIndices{};
    double maxScore = std::numeric_limits<double>::min();
    for (std::size_t ib = 0; ib < spacePointsOnTrack.size() - 2; ++ib) {
      for (std::size_t im = ib + 1; im < spacePointsOnTrack.size() - 1; ++im) {
        const double bmDeltaR =
            spacePointsOnTrack[im]->r() - spacePointsOnTrack[ib]->r();
        const double bmAbsDeltaZ =
            std::abs(spacePointsOnTrack[im]->z() - spacePointsOnTrack[ib]->z());
        if (bmDeltaR < 0) {
          ACTS_WARNING(
              "Space points are not sorted in r. Difference middle-bottom: "
              << bmDeltaR);
          continue;
        }
        if (bmDeltaR < m_cfg.deltaRMin || bmDeltaR > m_cfg.deltaRMax) {
          continue;
        }
        if (bmAbsDeltaZ < m_cfg.absDeltaZMin ||
            bmAbsDeltaZ > m_cfg.absDeltaZMax) {
          continue;
        }

        for (std::size_t it = im + 1; it < spacePointsOnTrack.size(); ++it) {
          const double mtDeltaR =
              spacePointsOnTrack[it]->r() - spacePointsOnTrack[im]->r();
          const double mtAbsDeltaZ = std::abs(spacePointsOnTrack[it]->z() -
                                              spacePointsOnTrack[im]->z());
          if (mtDeltaR < 0) {
            ACTS_WARNING(
                "Space points are not sorted in r. Difference top-middle: "
                << mtDeltaR);
            continue;
          }
          if (mtDeltaR < m_cfg.deltaRMin || mtDeltaR > m_cfg.deltaRMax) {
            continue;
          }
          if (mtAbsDeltaZ < m_cfg.absDeltaZMin ||
              mtAbsDeltaZ > m_cfg.absDeltaZMax) {
            continue;
          }

          const double score = bmDeltaR * mtDeltaR;

          if (score > maxScore) {
            seedFound = true;
            bestSPIndices = {ib, im, it};
            maxScore = score;
          }
        }
      }
    }

    if (seedFound) {
      SimSeed seed{*spacePointsOnTrack[bestSPIndices[0]],
                   *spacePointsOnTrack[bestSPIndices[1]],
                   *spacePointsOnTrack[bestSPIndices[2]]};

      Acts::ParticleHypothesis hypothesis =
          m_cfg.particleHypothesis.value_or(particle.hypothesis());

      seededParticles.insert(particle);
      tracks.emplace_back(std::move(track));
      seeds.emplace_back(seed);
      particleHypotheses.emplace_back(hypothesis);
    }
  }

  ACTS_VERBOSE("Found " << seeds.size() << " seeds");

  m_outputParticles(ctx, std::move(seededParticles));
  m_outputProtoTracks(ctx, std::move(tracks));
  m_outputSeeds(ctx, std::move(seeds));
  if (m_outputParticleHypotheses.isInitialized()) {
    m_outputParticleHypotheses(ctx, std::move(particleHypotheses));
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
