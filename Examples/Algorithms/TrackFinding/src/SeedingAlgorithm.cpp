// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"

#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/BinnedSPGroup.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/Seedfinder.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <stdexcept>

ActsExamples::SeedingAlgorithm::SeedingAlgorithm(
    ActsExamples::SeedingAlgorithm::Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("SeedingAlgorithm", lvl),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputSpacePoints.empty()) {
    throw std::invalid_argument("Missing space point input collections");
  }
  for (const auto& i : m_cfg.inputSpacePoints) {
    if (i.empty()) {
      throw std::invalid_argument("Invalid space point input collection");
    }
  }
  if (m_cfg.outputProtoTracks.empty()) {
    throw std::invalid_argument("Missing proto tracks output collection");
  }
  if (m_cfg.outputSeeds.empty()) {
    throw std::invalid_argument("Missing seeds output collection");
  }

  m_gridCfg.bFieldInZ = m_cfg.bFieldInZ;
  m_gridCfg.minPt = m_cfg.minPt;
  m_gridCfg.rMax = m_cfg.rMax;
  m_gridCfg.zMax = m_cfg.zMax;
  m_gridCfg.zMin = m_cfg.zMin;
  m_gridCfg.deltaRMax = m_cfg.deltaRMax;
  m_gridCfg.cotThetaMax = m_cfg.cotThetaMax;

  // construct seed filter
  Acts::SeedFilterConfig filterCfg;
  filterCfg.maxSeedsPerSpM = m_cfg.maxSeedsPerSpM;
  m_finderCfg.seedFilter = std::make_unique<Acts::SeedFilter<SimSpacePoint>>(
      Acts::SeedFilter<SimSpacePoint>(filterCfg));

  m_finderCfg.rMax = m_cfg.rMax;
  m_finderCfg.deltaRMin = m_cfg.deltaRMin;
  m_finderCfg.deltaRMax = m_cfg.deltaRMax;
  m_finderCfg.collisionRegionMin = m_cfg.collisionRegionMin;
  m_finderCfg.collisionRegionMax = m_cfg.collisionRegionMax;
  m_finderCfg.zMin = m_cfg.zMin;
  m_finderCfg.zMax = m_cfg.zMax;
  m_finderCfg.maxSeedsPerSpM = m_cfg.maxSeedsPerSpM;
  m_finderCfg.cotThetaMax = m_cfg.cotThetaMax;
  m_finderCfg.sigmaScattering = m_cfg.sigmaScattering;
  m_finderCfg.radLengthPerSeed = m_cfg.radLengthPerSeed;
  m_finderCfg.minPt = m_cfg.minPt;
  m_finderCfg.bFieldInZ = m_cfg.bFieldInZ;
  m_finderCfg.beamPos = Acts::Vector2(m_cfg.beamPosX, m_cfg.beamPosY);
  m_finderCfg.impactMax = m_cfg.impactMax;
}

ActsExamples::ProcessCode ActsExamples::SeedingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  // construct the combined input container of space point pointers from all
  // configured input sources.
  // pre-compute the total size required so we only need to allocate once
  size_t nSpacePoints = 0;
  for (const auto& isp : m_cfg.inputSpacePoints) {
    nSpacePoints += ctx.eventStore.get<SimSpacePointContainer>(isp).size();
  }
  std::vector<const SimSpacePoint*> spacePointPtrs;
  spacePointPtrs.reserve(nSpacePoints);
  for (const auto& isp : m_cfg.inputSpacePoints) {
    for (const auto& spacePoint :
         ctx.eventStore.get<SimSpacePointContainer>(isp)) {
      // since the event store owns the space points, their pointers should be
      // stable and we do not need to create local copies.
      spacePointPtrs.push_back(&spacePoint);
    }
  }

  // construct the seeding tools
  // covariance tool, extracts covariances per spacepoint as required
  auto extractGlobalQuantities =
      [=](const SimSpacePoint& sp, float, float,
          float) -> std::pair<Acts::Vector3, Acts::Vector2> {
    Acts::Vector3 position{sp.x(), sp.y(), sp.z()};
    Acts::Vector2 covariance{sp.varianceR(), sp.varianceZ()};
    return std::make_pair(position, covariance);
  };

  auto bottomBinFinder = std::make_shared<Acts::BinFinder<SimSpacePoint>>(
      Acts::BinFinder<SimSpacePoint>());
  auto topBinFinder = std::make_shared<Acts::BinFinder<SimSpacePoint>>(
      Acts::BinFinder<SimSpacePoint>());
  auto grid = Acts::SpacePointGridCreator::createGrid<SimSpacePoint>(m_gridCfg);
  auto spacePointsGrouping = Acts::BinnedSPGroup<SimSpacePoint>(
      spacePointPtrs.begin(), spacePointPtrs.end(), extractGlobalQuantities,
      bottomBinFinder, topBinFinder, std::move(grid), m_finderCfg);
  auto finder = Acts::Seedfinder<SimSpacePoint>(m_finderCfg);

  // run the seeding
  SimSeedContainer seeds;
  auto group = spacePointsGrouping.begin();
  auto groupEnd = spacePointsGrouping.end();
  for (; !(group == groupEnd); ++group) {
    const auto& groupSeeds =
        finder.createSeedsForGroup(group.bottom(), group.middle(), group.top());
    std::copy(groupSeeds.begin(), groupSeeds.end(), std::back_inserter(seeds));
  }

  // extract proto tracks, i.e. groups of measurement indices, from tracks seeds
  size_t nSeeds = seeds.size();
  ProtoTrackContainer protoTracks;
  protoTracks.reserve(nSeeds);
  for (const auto& seed : seeds) {
    ProtoTrack protoTrack;
    protoTrack.reserve(seed.sp().size());
    for (auto spacePointPtr : seed.sp()) {
      protoTrack.push_back(spacePointPtr->measurementIndex());
    }
    protoTracks.push_back(std::move(protoTrack));
  }

  ACTS_DEBUG("Created " << seeds.size() << " track seeds from "
                        << spacePointPtrs.size() << " space points");

  ctx.eventStore.add(m_cfg.outputSeeds, std::move(seeds));
  ctx.eventStore.add(m_cfg.outputProtoTracks, std::move(protoTracks));
  return ActsExamples::ProcessCode::SUCCESS;
}
