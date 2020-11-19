// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Seeding/SeedingAlgorithm.hpp"

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/BinnedSPGroup.hpp"
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/Seedfinder.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimVertex.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Io/Csv/CsvPlanarClusterReader.hpp"

#include <iostream>
#include <stdexcept>

using SimSpacePoint = ActsExamples::SimSpacePoint;

ActsExamples::SeedingAlgorithm::SeedingAlgorithm(
    ActsExamples::SeedingAlgorithm::Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("SeedingAlgorithm", lvl),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputClusters.empty()) {
    throw std::invalid_argument(
        "Missing clusters input collection with the hits");
  }
  if (m_cfg.outputSeeds.empty()) {
    throw std::invalid_argument("Missing output seeds collection");
  }
}

std::unique_ptr<SimSpacePoint> ActsExamples::SeedingAlgorithm::transformSP(
    std::size_t hit_id, const Acts::PlanarModuleCluster& cluster,
    const AlgorithmContext& ctx) const {
  const auto parameters = cluster.parameters();
  Acts::Vector2D localPos(parameters[0], parameters[1]);
  Acts::Vector3D globalPos(0, 0, 0);
  Acts::Vector3D globalFakeMom(1, 1, 1);

  // transform local into global position information
  globalPos = cluster.referenceObject().localToGlobal(ctx.geoContext, localPos,
                                                      globalFakeMom);
  auto cov = cluster.covariance();
  float x, y, z, r, varianceR, varianceZ;
  x = globalPos.x();
  y = globalPos.y();
  z = globalPos.z();
  r = std::sqrt(x * x + y * y);
  varianceR = cov(0, 0);
  varianceZ = cov(1, 1);
  std::unique_ptr<SimSpacePoint> sp(
      new SimSpacePoint{hit_id, x, y, z, r, varianceR, varianceZ});
  return sp;
}

ActsExamples::ProcessCode ActsExamples::SeedingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  Acts::SeedfinderConfig<SimSpacePoint> finderConf;
  // silicon detector max
  finderConf.rMax = m_cfg.rMax;
  finderConf.deltaRMin = m_cfg.deltaRMin;
  finderConf.deltaRMax = m_cfg.deltaRMax;
  finderConf.collisionRegionMin = m_cfg.collisionRegionMin;
  finderConf.collisionRegionMax = m_cfg.collisionRegionMax;
  finderConf.zMin = m_cfg.zMin;
  finderConf.zMax = m_cfg.zMax;
  finderConf.maxSeedsPerSpM = m_cfg.maxSeedsPerSpM;
  finderConf.cotThetaMax = m_cfg.cotThetaMax;
  finderConf.sigmaScattering = m_cfg.sigmaScattering;
  finderConf.radLengthPerSeed = m_cfg.radLengthPerSeed;
  finderConf.minPt = m_cfg.minPt;
  finderConf.bFieldInZ = m_cfg.bFieldInZ;
  finderConf.beamPos = m_cfg.beamPos;
  finderConf.impactMax = m_cfg.impactMax;

  // setup spacepoint grid config
  Acts::SpacePointGridConfig gridConf;
  gridConf.bFieldInZ = m_cfg.bFieldInZ;
  gridConf.minPt = m_cfg.minPt;
  gridConf.rMax = m_cfg.rMax;
  gridConf.zMax = m_cfg.zMax;
  gridConf.zMin = m_cfg.zMin;
  gridConf.deltaRMax = m_cfg.deltaRMax;
  gridConf.cotThetaMax = m_cfg.cotThetaMax;

  auto bottomBinFinder = std::make_shared<Acts::BinFinder<SimSpacePoint>>(
      Acts::BinFinder<SimSpacePoint>());
  auto topBinFinder = std::make_shared<Acts::BinFinder<SimSpacePoint>>(
      Acts::BinFinder<SimSpacePoint>());
  Acts::SeedFilterConfig sfconf;
  sfconf.maxSeedsPerSpM = m_cfg.maxSeedsPerSpM;
  finderConf.seedFilter = std::make_unique<Acts::SeedFilter<SimSpacePoint>>(
      Acts::SeedFilter<SimSpacePoint>(sfconf));
  Acts::Seedfinder<SimSpacePoint> seedFinder(finderConf);

  // covariance tool, sets covariances per spacepoint as required
  auto ct = [=](const SimSpacePoint& sp, float, float,
                float) -> Acts::Vector2D {
    return {sp.varianceR, sp.varianceZ};
  };

  const auto& clusters =
      ctx.eventStore
          .get<ActsExamples::GeometryIdMultimap<Acts::PlanarModuleCluster>>(
              m_cfg.inputClusters);

  // create the space points
  std::vector<const SimSpacePoint*> spVec;
  // since clusters are ordered, we simply count the hit_id as we read
  // clusters. Hit_id isn't stored in a cluster. This is how
  // CsvPlanarClusterWriter did it.
  std::size_t hit_id = 0;
  for (const auto& entry : clusters) {
    Acts::GeometryIdentifier geoId = entry.first;
    const Acts::PlanarModuleCluster& cluster = entry.second;
    std::size_t volumeId = geoId.volume();

    if (std::find(m_cfg.seedVolumes.begin(), m_cfg.seedVolumes.end(),
                  volumeId) != m_cfg.seedVolumes.end()) {
      auto sp = transformSP(hit_id, cluster, ctx).release();
      spVec.push_back(sp);
    }
    hit_id++;
  }

  // create grid with bin sizes according to the configured geometry
  std::unique_ptr<Acts::SpacePointGrid<SimSpacePoint>> grid =
      Acts::SpacePointGridCreator::createGrid<SimSpacePoint>(gridConf);
  auto spGroup = Acts::BinnedSPGroup<SimSpacePoint>(
      spVec.begin(), spVec.end(), ct, bottomBinFinder, topBinFinder,
      std::move(grid), finderConf);

  std::vector<std::vector<Acts::Seed<SimSpacePoint>>> seedVector;
  auto groupIt = spGroup.begin();
  auto endOfGroups = spGroup.end();
  for (; !(groupIt == endOfGroups); ++groupIt) {
    seedVector.push_back(seedFinder.createSeedsForGroup(
        groupIt.bottom(), groupIt.middle(), groupIt.top()));
  }

  int numSeeds = 0;
  for (auto& outVec : seedVector) {
    numSeeds += outVec.size();
  }

  ACTS_DEBUG(spVec.size() << " hits, " << seedVector.size() << " regions, "
                          << numSeeds << " seeds");

  ctx.eventStore.add(m_cfg.outputSeeds, std::move(seedVector));

  return ActsExamples::ProcessCode::SUCCESS;
}
