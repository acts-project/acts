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
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

ActsExamples::SeedingOrthogonalAlgorithm::SeedingOrthogonalAlgorithm(
    ActsExamples::SeedingOrthogonalAlgorithm::Config cfg,
    Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("SeedingAlgorithm", lvl),
      m_cfg(std::move(cfg)) {
  m_cfg.seedFilterConfig = m_cfg.seedFilterConfig.toInternalUnits();
  m_cfg.seedFinderConfig =
      m_cfg.seedFinderConfig.toInternalUnits().calculateDerivedQuantities();
  m_cfg.seedFinderOptions =
      m_cfg.seedFinderOptions.toInternalUnits().calculateDerivedQuantities(
          m_cfg.seedFinderConfig);

  printOptions();
  printConfig<SimSpacePoint>();

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
  if (m_cfg.seedFilterConfig.maxSeedsPerSpM !=
      m_cfg.seedFinderConfig.maxSeedsPerSpM) {
    throw std::invalid_argument("Inconsistent config maxSeedsPerSpM");
  }

  // construct seed filter
  m_cfg.seedFinderConfig.seedFilter =
      std::make_unique<Acts::SeedFilter<SimSpacePoint>>(
          Acts::SeedFilter<SimSpacePoint>(m_cfg.seedFilterConfig));

  m_finder = Acts::SeedFinderOrthogonal<SimSpacePoint>(m_cfg.seedFinderConfig);
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

  std::function<std::pair<Acts::Vector3, Acts::Vector2>(
      const SimSpacePoint *sp)>
      create_coordinates = [](const SimSpacePoint *sp) {
        Acts::Vector3 position(sp->x(), sp->y(), sp->z());
        Acts::Vector2 variance(sp->varianceR(), sp->varianceZ());
        return std::make_pair(position, variance);
      };

  SimSeedContainer seeds = finder.createSeeds(m_cfg.seedFinderOptions,
                                              spacePoints, create_coordinates);

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
      const IndexSourceLink &slink =
          spacePointPtr->sourceLinks()[0].get<IndexSourceLink>();
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

void ActsExamples::SeedingOrthogonalAlgorithm::printOptions() const {
  ACTS_DEBUG("SeedFinderOptions")
  ACTS_DEBUG("beamPos           " << m_cfg.seedFinderOptions.beamPos);
  // field induction
  ACTS_DEBUG("bFieldInZ         " << m_cfg.seedFinderOptions.bFieldInZ);
  // derived quantities
  ACTS_DEBUG("pTPerHelixRadius  " << m_cfg.seedFinderOptions.pTPerHelixRadius);
  ACTS_DEBUG("minHelixDiameter2 " << m_cfg.seedFinderOptions.minHelixDiameter2);
  ACTS_DEBUG("pT2perRadius      " << m_cfg.seedFinderOptions.pT2perRadius);
  ACTS_DEBUG("sigmapT2perRadius " << m_cfg.seedFinderOptions.sigmapT2perRadius);
  ACTS_DEBUG("...\n")
}

template <typename sp>
void ActsExamples::SeedingOrthogonalAlgorithm::printConfig() const {
  ACTS_DEBUG("SeedFinderOrthogonalConfig")
  ACTS_DEBUG("minPt                 " << m_cfg.seedFinderConfig.minPt);
  ACTS_DEBUG("deltaRMinTopSP        " << m_cfg.seedFinderConfig.deltaRMinTopSP);
  ACTS_DEBUG("deltaRMaxTopSP        " << m_cfg.seedFinderConfig.deltaRMaxTopSP);
  ACTS_DEBUG("deltaRMinBottomSP     "
             << m_cfg.seedFinderConfig.deltaRMinBottomSP);
  ACTS_DEBUG("deltaRMaxBottomSP     "
             << m_cfg.seedFinderConfig.deltaRMaxBottomSP);
  ACTS_DEBUG("impactMax             " << m_cfg.seedFinderConfig.impactMax);
  ACTS_DEBUG("maxPtScattering       "
             << m_cfg.seedFinderConfig.maxPtScattering);
  ACTS_DEBUG("collisionRegionMin    "
             << m_cfg.seedFinderConfig.collisionRegionMin);
  ACTS_DEBUG("collisionRegionMax    "
             << m_cfg.seedFinderConfig.collisionRegionMax);
  ACTS_DEBUG("zMin                  " << m_cfg.seedFinderConfig.zMin);
  ACTS_DEBUG("zMax                  " << m_cfg.seedFinderConfig.zMax);
  ACTS_DEBUG("rMax                  " << m_cfg.seedFinderConfig.rMax);
  ACTS_DEBUG("rMin                  " << m_cfg.seedFinderConfig.rMin);
  ACTS_DEBUG("highland              " << m_cfg.seedFinderConfig.highland);
  ACTS_DEBUG("maxScatteringAngle2   "
             << m_cfg.seedFinderConfig.maxScatteringAngle2);
  ACTS_DEBUG("...\n")
}
