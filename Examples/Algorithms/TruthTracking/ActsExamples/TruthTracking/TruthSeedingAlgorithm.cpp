// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/TruthSeedingAlgorithm.hpp"

#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <stdexcept>
#include <unordered_map>

ActsExamples::TruthSeedingAlgorithm::TruthSeedingAlgorithm(
    ActsExamples::TruthSeedingAlgorithm::Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("TruthSeedingAlgorithm", lvl),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input truth particles collection");
  }
  if (m_cfg.inputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument("Missing input hit-particles map collection");
  }
  if (m_cfg.inputSourceLinks.empty()) {
    throw std::invalid_argument("Missing source links input collection");
  }
  if (m_cfg.inputSpacePoints.empty()) {
    throw std::invalid_argument("Missing seeds or space point collection");
  }
  for (const auto& i : m_cfg.inputSpacePoints) {
    if (i.empty()) {
      throw std::invalid_argument("Invalid space point input collection");
    }
  }
  if (m_cfg.outputParticles.empty()) {
    throw std::invalid_argument("Missing output particles collection");
  }
  if (m_cfg.outputFullProtoTracks.empty()) {
    throw std::invalid_argument("Missing output full proto tracks collection");
  }
  if (m_cfg.outputSeeds.empty()) {
    throw std::invalid_argument("Missing seeds output collections");
  }
  if (m_cfg.outputProtoTracks.empty()) {
    throw std::invalid_argument("Missing proto tracks output collections");
  }
}

ActsExamples::ProcessCode ActsExamples::TruthSeedingAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  using HitParticlesMap = IndexMultimap<ActsFatras::Barcode>;

  // prepare input collections
  const auto& particles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);
  const auto& hitParticlesMap =
      ctx.eventStore.get<HitParticlesMap>(m_cfg.inputMeasurementParticlesMap);
  // compute particle_id -> {hit_id...} map from the
  // hit_id -> {particle_id...} map on the fly.
  const auto& particleHitsMap = invertIndexMultimap(hitParticlesMap);

  SimSpacePointContainer spacePoints;
  for (const auto& isp : m_cfg.inputSpacePoints) {
    const auto& sps = ctx.eventStore.get<SimSpacePointContainer>(isp);
    std::copy(sps.begin(), sps.end(), std::back_inserter(spacePoints));
  }

  SimParticleContainer seededParticles;
  ProtoTrackContainer fullTracks;
  SimSeedContainer seeds;
  ProtoTrackContainer tracks;

  seededParticles.reserve(particles.size());
  fullTracks.reserve(particles.size());
  seeds.reserve(particles.size());
  tracks.reserve(particles.size());

  std::unordered_map<Index, const SimSpacePoint*> spMap;

  for (const SimSpacePoint& sp : spacePoints) {
    if (sp.sourceLinks().empty()) {
      ACTS_WARNING("Missing source link in space point");
      continue;
    }
    for (const auto& slink : sp.sourceLinks()) {
      const IndexSourceLink& islink = slink.get<IndexSourceLink>();
      spMap.emplace(islink.index(), &sp);
    }
  }

  for (const auto& particle : particles) {
    // find the corresponding hits for this particle
    const auto& hits =
        makeRange(particleHitsMap.equal_range(particle.particleId()));
    // fill hit indices to create the proto track
    ProtoTrack fullTrack;
    fullTrack.reserve(hits.size());
    for (const auto& hit : hits) {
      fullTrack.emplace_back(hit.second);
    }

    // The list of hits and the initial start parameters
    if (fullTrack.size() < 3) {
      ACTS_WARNING("Particle " << particle << " has less than 3 hits");
      continue;
    }
    // Space points on the proto track
    std::vector<const SimSpacePoint*> spacePointsOnTrack;
    spacePointsOnTrack.reserve(fullTrack.size());
    // Loop over the hit index on the proto track to find the space points
    for (const auto& hitIndex : fullTrack) {
      auto it = spMap.find(hitIndex);
      if (it != spMap.end()) {
        spacePointsOnTrack.push_back(it->second);
      }
    }
    // At least three space points are required
    if (spacePointsOnTrack.size() < 3) {
      continue;
    }
    // Sort the space points
    std::sort(spacePointsOnTrack.begin(), spacePointsOnTrack.end(),
              [](const SimSpacePoint* lhs, const SimSpacePoint* rhs) {
                return std::hypot(lhs->r(), lhs->z()) <
                       std::hypot(rhs->r(), rhs->z());
              });

    // Loop over the found space points to find the seed with maxium deltaR
    // betweent the bottom and top space point
    // @todo add the check of deltaZ
    bool seedFound = false;
    std::array<size_t, 3> bestSPIndices{};
    double maxDeltaR = std::numeric_limits<double>::min();
    for (size_t ib = 0; ib < spacePointsOnTrack.size() - 2; ++ib) {
      for (size_t im = ib + 1; im < spacePointsOnTrack.size() - 1; ++im) {
        for (size_t it = im + 1; it < spacePointsOnTrack.size(); ++it) {
          double bmDeltaR = std::abs(spacePointsOnTrack[im]->r() -
                                     spacePointsOnTrack[ib]->r());
          double mtDeltaR = std::abs(spacePointsOnTrack[it]->r() -
                                     spacePointsOnTrack[im]->r());
          if (bmDeltaR >= m_cfg.deltaRMin and bmDeltaR <= m_cfg.deltaRMax and
              mtDeltaR >= m_cfg.deltaRMin and mtDeltaR <= m_cfg.deltaRMax) {
            if ((bmDeltaR + mtDeltaR) > maxDeltaR) {
              maxDeltaR = bmDeltaR + mtDeltaR;
              bestSPIndices = {ib, im, it};
              seedFound = true;
            }
          }
        }
      }
    }

    if (seedFound) {
      SimSeed seed{
          *spacePointsOnTrack[bestSPIndices[0]],
          *spacePointsOnTrack[bestSPIndices[1]],
          *spacePointsOnTrack[bestSPIndices[2]],
          static_cast<float>(spacePointsOnTrack[bestSPIndices[1]]->z())};

      ProtoTrack track;
      track.reserve(3);
      for (const auto& sp : seed.sp()) {
        if (sp->sourceLinks().empty()) {
          ACTS_ERROR("Missing source link in the space point")
          continue;
        }
        const auto slink = sp->sourceLinks()[0].get<IndexSourceLink>();
        track.push_back(slink.index());
      }

      seededParticles.insert(particle);
      fullTracks.emplace_back(std::move(fullTrack));
      seeds.push_back(seed);
      tracks.emplace_back(track);
    }
  }

  ACTS_VERBOSE("Found " << seeds.size() << " seeds");

  ctx.eventStore.add(m_cfg.outputParticles, std::move(seededParticles));
  ctx.eventStore.add(m_cfg.outputFullProtoTracks, std::move(fullTracks));
  ctx.eventStore.add(m_cfg.outputSeeds, std::move(seeds));
  ctx.eventStore.add(m_cfg.outputProtoTracks, std::move(tracks));
  return ActsExamples::ProcessCode::SUCCESS;
}
