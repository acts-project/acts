// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Seeding/SeedingAlgorithm.hpp"

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/BinnedSPGroup.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/Seedfinder.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimVertex.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <stdexcept>

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
  if (not m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
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

std::unique_ptr<ActsExamples::SimSpacePoint>
ActsExamples::SeedingAlgorithm::makeSpacePoint(
    const Measurement& measurement, const Acts::Surface& surface,
    const Acts::GeometryContext& geoCtx) const {
  SimSpacePoint sp = std::visit(
      [&](const auto& m) {
        // since we do not know if and where the local parameters are contained
        // in the measurement, we always transform to the bound space where we
        // do not their location. if the local parameters are not measured, this
        // results in a zero location.
        Acts::BoundVector localPar = m.expander() * m.parameters();
        Acts::BoundSymMatrix localCov =
            m.expander() * m.covariance() * m.expander().transpose();

        // transform local position to global coordinates
        Acts::Vector3D globalFakeMom(1, 1, 1);
        Acts::Vector3D globalPos = surface.localToGlobal(
            geoCtx, {localPar[Acts::eBoundLoc0], localPar[Acts::eBoundLoc1]},
            globalFakeMom);

        // transfrom local covariance to global coordinates
        Acts::BoundToFreeMatrix jacToGlobal =
            surface.jacobianLocalToGlobal(geoCtx, localPar);
        Acts::SymMatrix3D globalCov =
            (jacToGlobal * localCov * jacToGlobal.transpose())
                .block<3, 3>(Acts::eFreePos0, Acts::eFreePos0);

        // create space point
        float x = globalPos[Acts::ePos0];
        float y = globalPos[Acts::ePos1];
        float z = globalPos[Acts::ePos2];
        float r = std::hypot(x, y);
        float varianceR = globalCov(0, 0);
        float varianceZ = globalCov(1, 1);
        return SimSpacePoint{
            m.sourceLink().index(), x, y, z, r, varianceR, varianceZ};
      },
      measurement);
  return std::make_unique<SimSpacePoint>(std::move(sp));
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
  for (Index imeas = 0u; imeas < measurements.size(); ++imeas) {
    const auto& meas = measurements[imeas];
    const auto geoId = std::visit(
        [](const auto& m) { return m.sourceLink().geometryId(); }, meas);
    const auto volumeId = geoId.volume();
    const auto layerId = geoId.layer();

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
    } else {
      continue;
    }

    // lookup surface
    const Acts::Surface* surface = m_cfg.trackingGeometry->findSurface(geoId);
    if (not surface) {
      ACTS_ERROR("Could not find surface " << geoId);
      return ProcessCode::ABORT;
    }

    auto sp = makeSpacePoint(meas, *surface, ctx.geoContext).release();
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
