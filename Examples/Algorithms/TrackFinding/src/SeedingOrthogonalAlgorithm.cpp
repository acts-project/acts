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

void printOptions( const Acts::SeedFinderOptions& options) {
  std::cout << "SeedFinderOptions\n"; 
  std::cout << "beamPos           "  << options.beamPos << std::endl;
  // field induction
  std::cout << "bFieldInZ         "  << options.bFieldInZ << std::endl;
  // derived quantities
  std::cout << "pTPerHelixRadius  "  << options.pTPerHelixRadius << std::endl;
  std::cout << "minHelixDiameter2 "  << options.minHelixDiameter2 << std::endl;
  std::cout << "pT2perRadius      "  << options.pT2perRadius << std::endl;
  std::cout << "sigmapT2perRadius "  << options.sigmapT2perRadius << std::endl;
  std::cout << "...\n"; 

}

template<typename sp>
void printConfig( const Acts::SeedFinderOrthogonalConfig<sp>& config ) {
  std::cout << "SeedFinderOrthogonalConfig\n"; 

    std::cout << "minPt                 " << config.minPt  << std::endl;
    std::cout << "deltaRMinTopSP        " << config.deltaRMinTopSP << std::endl;
    std::cout << "deltaRMaxTopSP        " << config.deltaRMaxTopSP << std::endl;
    std::cout << "deltaRMinBottomSP     " << config.deltaRMinBottomSP << std::endl;
    std::cout << "deltaRMaxBottomSP     " << config.deltaRMaxBottomSP << std::endl;
    std::cout << "impactMax             " << config.impactMax << std::endl;
    std::cout << "maxPtScattering       " << config.maxPtScattering  << std::endl;
    std::cout << "collisionRegionMin    " << config.collisionRegionMin << std::endl;
    std::cout << "collisionRegionMax    " << config.collisionRegionMax << std::endl;
    std::cout << "zMin                  " << config.zMin << std::endl;
    std::cout << "zMax                  " << config.zMax << std::endl;
    std::cout << "rMax                  " << config.rMax << std::endl;
    std::cout << "rMin                  " << config.rMin << std::endl;

    std::cout << "highland              " << config.highland << std::endl;
    std::cout << "maxScatteringAngle2   " << config.maxScatteringAngle2 << std::endl;
  std::cout << "...\n"; 

}

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

  printOptions(m_cfg.seedFinderOptions);
  printConfig(m_cfg.seedFinderConfig);

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

  SimSeedContainer seeds =
      m_finder.createSeeds(m_cfg.seedFinderOptions, spacePoints);

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
