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
#include "ActsExamples/EventData/SimSeed.hpp"
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
  if (m_cfg.outputTrackParametersSeedMap.empty()) {
    throw std::invalid_argument(
        "Missing output trackparameters-to-seed collection");
  }
  if (not m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }

  // Set up the track parameters covariance (the same for all tracks)
  m_covariance(Acts::eBoundLoc0, Acts::eBoundLoc0) =
      cfg.sigmaLoc0 * m_cfg.sigmaLoc0;
  m_covariance(Acts::eBoundLoc1, Acts::eBoundLoc1) =
      cfg.sigmaLoc1 * m_cfg.sigmaLoc1;
  m_covariance(Acts::eBoundPhi, Acts::eBoundPhi) =
      cfg.sigmaPhi * m_cfg.sigmaPhi;
  m_covariance(Acts::eBoundTheta, Acts::eBoundTheta) =
      cfg.sigmaTheta * m_cfg.sigmaTheta;
  m_covariance(Acts::eBoundQOverP, Acts::eBoundQOverP) =
      cfg.sigmaQOverP * m_cfg.sigmaQOverP;
  m_covariance(Acts::eBoundTime, Acts::eBoundTime) =
      m_cfg.sigmaT0 * m_cfg.sigmaT0;
}

Acts::Vector3 ActsExamples::TrackParamsEstimationAlgorithm::getField(
    const Acts::Vector3& position) const {
  return std::visit(
      [&](const auto& inputField) -> Acts::Vector3 {
        return inputField->getField(position);
      },
      m_cfg.magneticField);
}

ActsExamples::ProcessCode ActsExamples::TrackParamsEstimationAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // @todo storing the seed as a ProtoTrack in the SeedingAlgorithm. The input
  // will be ProtoTrackContainer then.
  const auto& seeds = ctx.eventStore.get<SimSeedContainer>(m_cfg.inputSeeds);
  // need source links to get the geometry identifer
  const auto& sourceLinks =
      ctx.eventStore.get<IndexSourceLinkContainer>(m_cfg.inputSourceLinks);

  TrackParametersContainer trackParameters;
  TrackParametersSeedMap trackParametersSeedMap;
  trackParameters.reserve(seeds.size());

  for (Index iregion = 0; iregion < seeds.size(); ++iregion) {
    const auto& regionSeeds = seeds[iregion];
    for (Index iseed = 0; iseed < regionSeeds.size(); ++iseed) {
      const auto& seed = regionSeeds[iseed];
      // Get the bottom space point and its reference surface
      // @todo do we need to sort the sps first
      const auto bottomSP = seed.sp().front();
      const auto hitIdx = bottomSP->measurementIndex();
      const auto sourceLink = sourceLinks.nth(hitIdx);
      auto geoId = sourceLink->geometryId();
      const Acts::Surface* surface = m_cfg.trackingGeometry->findSurface(geoId);
      if (surface == nullptr) {
        ACTS_WARNING("surface with geoID "
                     << geoId << " is not found in the tracking gemetry");
        continue;
      }

      // Get the magnetic field at the bottom space point
      Acts::Vector3 field = m_cfg.magneticField->getField(
          {bottomSP->x(), bottomSP->y(), bottomSP->z()});

      auto optParams = Acts::estimateTrackParamsFromSeed(
          ctx.geoContext, seed.sp().begin(), seed.sp().end(), *surface, field,
          m_cfg.bFieldMin);
      if (not optParams.has_value()) {
        ACTS_WARNING("Estimation of track parameters from seed "
                     << iseed << " in region " << iregion << " failed.");
        continue;
      } else {
        const auto& params = optParams.value();
        double charge = std::copysign(1, params[Acts::eBoundQOverP]);
        trackParameters.emplace_back(surface->getSharedPtr(), params, charge,
                                     m_covariance);
        trackParametersSeedMap.emplace(trackParameters.size() - 1,
                                       GroupedSeedIndex{iregion, iseed});
      }
    }
  }

  ctx.eventStore.add(m_cfg.outputTrackParameters, std::move(trackParameters));
  ctx.eventStore.add(m_cfg.outputTrackParametersSeedMap,
                     std::move(trackParametersSeedMap));
  return ActsExamples::ProcessCode::SUCCESS;
}
