// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/BinnedSPGroup.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/Seedfinder.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
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
  m_finderCfg.beamPos = m_cfg.beamPos;
  m_finderCfg.impactMax = m_cfg.impactMax;
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
        Acts::Vector3 globalFakeMom(1, 1, 1);
        Acts::Vector3 globalPos = surface.localToGlobal(
            geoCtx, {localPar[Acts::eBoundLoc0], localPar[Acts::eBoundLoc1]},
            globalFakeMom);

        // transfrom local covariance to global coordinates
        Acts::BoundToFreeMatrix jacToGlobal =
            surface.jacobianLocalToGlobal(geoCtx, localPar);
        Acts::SymMatrix3 globalCov =
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

  // covariance tool, sets covariances per spacepoint as required
  auto ct = [=](const SimSpacePoint& sp, float, float, float) -> Acts::Vector2 {
    return {sp.varianceR(), sp.varianceZ()};
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
  auto grid = Acts::SpacePointGridCreator::createGrid<SimSpacePoint>(m_gridCfg);
  Acts::BinnedSPGroup<SimSpacePoint> spGroup(spVec.begin(), spVec.end(), ct,
                                             bottomBinFinder, topBinFinder,
                                             std::move(grid), m_finderCfg);
  Acts::Seedfinder<SimSpacePoint> finder(m_finderCfg);

  std::vector<std::vector<Acts::Seed<SimSpacePoint>>> seedVector;
  auto groupIt = spGroup.begin();
  auto endOfGroups = spGroup.end();
  for (; !(groupIt == endOfGroups); ++groupIt) {
    seedVector.push_back(finder.createSeedsForGroup(
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
