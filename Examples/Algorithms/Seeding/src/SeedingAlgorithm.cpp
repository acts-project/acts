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

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimVertex.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include <iostream>
#include <stdexcept>

using SimSpacePoint = ActsExamples::SimSpacePoint;
using ConcreteMeasurement =
    Acts::Measurement<ActsExamples::IndexSourceLink, Acts::BoundIndices,
                      Acts::eBoundLoc0, Acts::eBoundLoc1>;

ActsExamples::SeedingAlgorithm::SeedingAlgorithm(
    ActsExamples::SeedingAlgorithm::Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("SeedingAlgorithm", lvl),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument(
        "Missing measurement input collection with the hits");
  }
  if (m_cfg.outputSeeds.empty()) {
    throw std::invalid_argument("Missing output seeds collection");
  }

  // silicon detector max
  m_cfg.finderConf.rMax = m_cfg.rMax;
  m_cfg.finderConf.deltaRMin = m_cfg.deltaRMin;
  m_cfg.finderConf.deltaRMax = m_cfg.deltaRMax;
  m_cfg.finderConf.collisionRegionMin = m_cfg.collisionRegionMin;
  m_cfg.finderConf.collisionRegionMax = m_cfg.collisionRegionMax;
  m_cfg.finderConf.zMin = m_cfg.zMin;
  m_cfg.finderConf.zMax = m_cfg.zMax;
  m_cfg.finderConf.maxSeedsPerSpM = m_cfg.maxSeedsPerSpM;
  m_cfg.finderConf.cotThetaMax = m_cfg.cotThetaMax;
  m_cfg.finderConf.sigmaScattering = m_cfg.sigmaScattering;
  m_cfg.finderConf.radLengthPerSeed = m_cfg.radLengthPerSeed;
  m_cfg.finderConf.minPt = m_cfg.minPt;
  m_cfg.finderConf.bFieldInZ = m_cfg.bFieldInZ;
  m_cfg.finderConf.beamPos = m_cfg.beamPos;
  m_cfg.finderConf.impactMax = m_cfg.impactMax;

  m_cfg.sfconf.maxSeedsPerSpM = m_cfg.maxSeedsPerSpM;
  m_cfg.finderConf.seedFilter =
      std::make_unique<Acts::SeedFilter<SimSpacePoint>>(
          Acts::SeedFilter<SimSpacePoint>(m_cfg.sfconf));
  m_cfg.gridConf.bFieldInZ = m_cfg.bFieldInZ;
  m_cfg.gridConf.minPt = m_cfg.minPt;
  m_cfg.gridConf.rMax = m_cfg.rMax;
  m_cfg.gridConf.zMax = m_cfg.zMax;
  m_cfg.gridConf.zMin = m_cfg.zMin;
  m_cfg.gridConf.deltaRMax = m_cfg.deltaRMax;
  m_cfg.gridConf.cotThetaMax = m_cfg.cotThetaMax;
}

std::unique_ptr<SimSpacePoint> ActsExamples::SeedingAlgorithm::transformSP(
    const unsigned int hit_id, const ConcreteMeasurement meas,
    const AlgorithmContext& ctx) const {
  const auto parameters = meas.parameters();
  Acts::Vector2D localPos(parameters[0], parameters[1]);
  Acts::Vector3D globalPos(0, 0, 0);
  Acts::Vector3D globalFakeMom(1, 1, 1);

  // transform local into global position information
  const auto& surf = meas.referenceObject();

  globalPos = meas.referenceObject().localToGlobal(ctx.geoContext, localPos,
                                                   globalFakeMom);
  auto cov_local = meas.covariance();
  auto rframe = surf.referenceFrame(ctx.geoContext, globalPos, globalFakeMom);
  auto jacToGlobal = rframe.topLeftCorner<3, 2>();
  auto cov_global = jacToGlobal * cov_local * jacToGlobal.transpose();

  float x, y, z, r, varianceR, varianceZ;
  x = globalPos.x();
  y = globalPos.y();
  z = globalPos.z();
  r = std::sqrt(x * x + y * y);
  varianceR = cov_global(0, 0);
  varianceZ = cov_global(1, 1);
  std::unique_ptr<SimSpacePoint> sp(
      new SimSpacePoint{hit_id, x, y, z, r, varianceR, varianceZ});

  return sp;
}

ActsExamples::ProcessCode ActsExamples::SeedingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  auto bottomBinFinder = std::make_shared<Acts::BinFinder<SimSpacePoint>>(
      Acts::BinFinder<SimSpacePoint>());
  auto topBinFinder = std::make_shared<Acts::BinFinder<SimSpacePoint>>(
      Acts::BinFinder<SimSpacePoint>());

  Acts::Seedfinder<SimSpacePoint> seedFinder(m_cfg.finderConf);

  // covariance tool, sets covariances per spacepoint as required
  auto ct = [=](const SimSpacePoint& sp, float, float,
                float) -> Acts::Vector2D {
    return {sp.varianceR, sp.varianceZ};
  };

  const auto& measurements =
      ctx.eventStore.get<MeasurementContainer>(m_cfg.inputMeasurements);

  // create the space points
  std::vector<const SimSpacePoint*> spVec;
  unsigned int hit_id = 0;
  for (const auto& fullMeas : measurements) {
    const auto& meas = std::get<ConcreteMeasurement>(fullMeas);
    const auto& surface = meas.referenceObject();
    const auto& geoId = surface.geometryId();
    unsigned int volumeId = geoId.volume();
    unsigned int layerId = geoId.layer();

    // volumes and layers for seed finding
    if (volumeId == m_cfg.barrelVolume) {
      if (std::find(m_cfg.barrelLayers.begin(), m_cfg.barrelLayers.end(),
                    layerId) == m_cfg.barrelLayers.end())
        continue;
    } else if (volumeId == m_cfg.posEndcapVolume) {
      if (std::find(m_cfg.posEndcapLayers.begin(), m_cfg.posEndcapLayers.end(),
                    layerId) == m_cfg.posEndcapLayers.end())
        continue;
    } else if (volumeId == m_cfg.negEndcapVolume) {
      if (std::find(m_cfg.negEndcapLayers.begin(), m_cfg.negEndcapLayers.end(),
                    layerId) == m_cfg.negEndcapLayers.end())
        continue;
    } else
      continue;

    auto sp = transformSP(hit_id, meas, ctx).release();
    spVec.push_back(sp);
    hit_id++;
  }
  Acts::SeedfinderConfig<SimSpacePoint> fconf = m_cfg.finderConf;
  std::unique_ptr<Acts::SpacePointGrid<SimSpacePoint>> grid =
      Acts::SpacePointGridCreator::createGrid<SimSpacePoint>(m_cfg.gridConf);

  auto spGroup = Acts::BinnedSPGroup<SimSpacePoint>(
      spVec.begin(), spVec.end(), ct, bottomBinFinder, topBinFinder,
      std::move(grid), fconf);

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
