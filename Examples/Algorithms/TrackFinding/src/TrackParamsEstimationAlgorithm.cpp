// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Seeding/EstimateTrackParamsFromSeed.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <map>
#include <stdexcept>

ActsExamples::TrackParamsEstimationAlgorithm::TrackParamsEstimationAlgorithm(
    ActsExamples::TrackParamsEstimationAlgorithm::Config cfg,
    Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("TrackParamsEstimationAlgorithm", lvl),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputSeeds.empty()) {
    throw std::invalid_argument("Missing input seeds collection");
  }
  if (m_cfg.inputSourceLinks.empty()) {
    throw std::invalid_argument("Missing input source links collection");
  }
  if (m_cfg.outputTrackParameters.empty()) {
    throw std::invalid_argument("Missing output track parameters collection");
  }
  if (m_cfg.outputTrackParamsSeedMap.empty()) {
    throw std::invalid_argument(
        "Missing output trackparameters-to-seed collection");
  }
  if (not m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }
}

ActsExamples::ProcessCode ActsExamples::TrackParamsEstimationAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  using SeedContainer = std::vector<std::vector<Acts::Seed<SimSpacePoint>>>;
  using TrackParamsSeedMap = std::map<Index, std::pair<Index, Index>>;
  const auto& seeds = ctx.eventStore.get<SeedContainer>(m_cfg.inputSeeds);
  // need source links to get the geometry identifer
  const auto& sourceLinks =
      ctx.eventStore.get<IndexSourceLinkContainer>(m_cfg.inputSourceLinks);

  TrackParametersContainer trackParameters;
  TrackParamsSeedMap trackParamsSeedMap;
  trackParameters.reserve(seeds.size());

  for (size_t iregion = 0; iregion < seeds.size(); ++iregion) {
    const auto& regionSeeds = seeds[iregion];
    for (size_t iseed = 0; iseed < regionSeeds.size(); ++iseed) {
      const auto& seed = regionSeeds[iseed];
      // Get the bottom space point and its reference surface
      // @todo do we need to sort the sps first
      const auto bottomSP = seed.sp().front();
      const auto hitIdx = bottomSP->measurementIndex();
      const auto sourceLink = sourceLinks.nth(hitIdx);
      auto geoId = sourceLink->geometryId();
      const Acts::Surface* surface = m_cfg.trackingGeometry->findSurface(geoId);
      if (not surface) {
        ACTS_WARNING("surface with geoID "
                     << geoId << " is not found in the tracking gemetry");
        continue;
      }

      // Get the magnetic field at the bottom space point
      Acts::Vector3 field = m_cfg.bFieldGetter(
          Acts::Vector3(bottomSP->x(), bottomSP->y(), bottomSP->z()));
      // Estimate the track parameters from seed
      auto optParams = Acts::estimateTrackParamsFromSeed(
          ctx.geoContext, seed.sp().begin(), seed.sp().end(), *surface, field,
          m_cfg.bFieldMin);
      if (not optParams.has_value()) {
        ACTS_WARNING("Estimation of track parameters from seed "
                     << iseed << " in region " << iregion << " failed.");
        continue;
      } else {
        const auto& params = optParams.value();
        double charge =
            std::abs(params[Acts::eBoundQOverP]) / params[Acts::eBoundQOverP];

        // compute momentum-dependent resolutions
        const double& sigmaLoc0 = m_cfg.sigmaLoc0;
        const double& sigmaLoc1 = m_cfg.sigmaLoc1;
        const double& sigmaPhi = m_cfg.sigmaPhi;
        const double& sigmaTheta = m_cfg.sigmaTheta;
        const double& sigmaQOverP = m_cfg.sigmaQOverP;
        const double& sigmaT0 = m_cfg.sigmaT0;

        Acts::BoundSymMatrix cov = Acts::BoundSymMatrix::Zero();
        cov(Acts::eBoundLoc0, Acts::eBoundLoc0) = sigmaLoc0 * sigmaLoc0;
        cov(Acts::eBoundLoc1, Acts::eBoundLoc1) = sigmaLoc1 * sigmaLoc1;
        cov(Acts::eBoundPhi, Acts::eBoundPhi) = sigmaPhi * sigmaPhi;
        cov(Acts::eBoundTheta, Acts::eBoundTheta) = sigmaTheta * sigmaTheta;
        cov(Acts::eBoundQOverP, Acts::eBoundQOverP) = sigmaQOverP * sigmaQOverP;
        cov(Acts::eBoundTime, Acts::eBoundTime) = sigmaT0 * sigmaT0;

        trackParameters.emplace_back(surface->getSharedPtr(), params, charge,
                                     cov);
        trackParamsSeedMap.emplace(trackParameters.size() - 1,
                                   std::make_pair(iregion, iseed));
      }
    }
  }

  ctx.eventStore.add(m_cfg.outputTrackParameters, std::move(trackParameters));
  ctx.eventStore.add(m_cfg.outputTrackParamsSeedMap,
                     std::move(trackParamsSeedMap));
  return ActsExamples::ProcessCode::SUCCESS;
}
