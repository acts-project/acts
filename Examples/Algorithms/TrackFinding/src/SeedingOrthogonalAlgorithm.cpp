// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/SeedingOrthogonalAlgorithm.hpp"

#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/SeedFinderOrthogonal.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

ActsExamples::SeedingOrthogonalAlgorithm::SeedingOrthogonalAlgorithm(
    ActsExamples::SeedingOrthogonalAlgorithm::Config cfg,
    Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("SeedingAlgorithm", lvl),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputSpacePoints.empty()) {
    throw std::invalid_argument("Missing space point input collections");
  }
  for (const auto &i : m_cfg.inputSpacePoints) {
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

  // construct seed filter
  Acts::SeedFilterConfig filterCfg;
  filterCfg.maxSeedsPerSpM = m_cfg.maxSeedsPerSpM;
  m_cfg.seedFinderConfig.seedFilter =
      std::make_unique<Acts::SeedFilter<SimSpacePoint>>(
          Acts::SeedFilter<SimSpacePoint>(filterCfg));

  m_cfg.seedFinderConfig.rMax = m_cfg.rMax;
  m_cfg.seedFinderConfig.deltaRMinTopSP = m_cfg.deltaRMinTopSP;
  m_cfg.seedFinderConfig.deltaRMaxTopSP = m_cfg.deltaRMaxTopSP;
  m_cfg.seedFinderConfig.deltaRMinBottomSP = m_cfg.deltaRMinBottomSP;
  m_cfg.seedFinderConfig.deltaRMaxBottomSP = m_cfg.deltaRMaxBottomSP;
  m_cfg.seedFinderConfig.collisionRegionMin = m_cfg.collisionRegionMin;
  m_cfg.seedFinderConfig.collisionRegionMax = m_cfg.collisionRegionMax;
  m_cfg.seedFinderConfig.zMin = m_cfg.zMin;
  m_cfg.seedFinderConfig.zMax = m_cfg.zMax;
  m_cfg.seedFinderConfig.maxSeedsPerSpM = m_cfg.maxSeedsPerSpM;
  m_cfg.seedFinderConfig.cotThetaMax = m_cfg.cotThetaMax;
  m_cfg.seedFinderConfig.sigmaScattering = m_cfg.sigmaScattering;
  m_cfg.seedFinderConfig.radLengthPerSeed = m_cfg.radLengthPerSeed;
  m_cfg.seedFinderConfig.minPt = m_cfg.minPt;
  m_cfg.seedFinderConfig.bFieldInZ = m_cfg.bFieldInZ;
  m_cfg.seedFinderConfig.beamPos =
      Acts::Vector2(m_cfg.beamPosX, m_cfg.beamPosY);
  m_cfg.seedFinderConfig.impactMax = m_cfg.impactMax;

  // calculation of scattering using the highland formula
  // convert pT to p once theta angle is known
  m_cfg.seedFinderConfig.highland =
      13.6 * std::sqrt(m_cfg.seedFinderConfig.radLengthPerSeed) *
      (1 + 0.038 * std::log(m_cfg.seedFinderConfig.radLengthPerSeed));
  float maxScatteringAngle =
      m_cfg.seedFinderConfig.highland / m_cfg.seedFinderConfig.minPt;
  m_cfg.seedFinderConfig.maxScatteringAngle2 =
      maxScatteringAngle * maxScatteringAngle;
  // helix radius in homogeneous magnetic field. Units are Kilotesla, MeV and
  // millimeter
  // TODO: change using ACTS units
  m_cfg.seedFinderConfig.pTPerHelixRadius =
      300. * m_cfg.seedFinderConfig.bFieldInZ;
  m_cfg.seedFinderConfig.minHelixDiameter2 =
      std::pow(m_cfg.seedFinderConfig.minPt * 2 /
                   m_cfg.seedFinderConfig.pTPerHelixRadius,
               2);

  m_cfg.seedFinderConfig.pT2perRadius = std::pow(
      m_cfg.seedFinderConfig.highland / m_cfg.seedFinderConfig.pTPerHelixRadius,
      2);
}

ActsExamples::ProcessCode ActsExamples::SeedingOrthogonalAlgorithm::execute(
    const AlgorithmContext &ctx) const {
  std::vector<const SimSpacePoint *> spacePoints;

  for (const auto &isp : m_cfg.inputSpacePoints) {
    for (const auto &spacePoint :
         ctx.eventStore.get<SimSpacePointContainer>(isp)) {
      spacePoints.push_back(&spacePoint);
    }
  }

  Acts::SeedFinderOrthogonal<SimSpacePoint> finder(m_cfg.seedFinderConfig);

  SimSeedContainer seeds = finder.createSeeds(spacePoints);

  // extract proto tracks, i.e. groups of measurement indices, from tracks seeds
  size_t nSeeds = seeds.size();
  ProtoTrackContainer protoTracks;
  protoTracks.reserve(nSeeds);
  for (const auto &seed : seeds) {
    ProtoTrack protoTrack;
    protoTrack.reserve(seed.sp().size());
    for (auto spacePointPtr : seed.sp()) {
      if (spacePointPtr->sourceLinks().empty()) {
        ACTS_WARNING("Missing sourcelink in space point");
        continue;
      }
      const auto slink = static_cast<const IndexSourceLink &>(
          *(spacePointPtr->sourceLinks()[0]));
      protoTrack.push_back(slink.index());
    }
    protoTracks.push_back(std::move(protoTrack));
  }

  ACTS_DEBUG("Created " << seeds.size() << " track seeds from "
                        << spacePoints.size() << " space points");

  ctx.eventStore.add(m_cfg.outputSeeds, std::move(seeds));
  ctx.eventStore.add(m_cfg.outputProtoTracks, std::move(protoTracks));

  return ActsExamples::ProcessCode::SUCCESS;
}
