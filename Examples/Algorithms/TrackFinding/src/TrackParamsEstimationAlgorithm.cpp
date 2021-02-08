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
#include "ActsExamples/EventData/ProtoTrack.hpp"
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
  if (m_cfg.inputSpacePoints.empty()) {
    throw std::invalid_argument("Missing space point input collections");
  }
  for (const auto& i : m_cfg.inputSpacePoints) {
    if (i.empty()) {
      throw std::invalid_argument("Invalid space point input collection");
    }
  }
  if (m_cfg.inputProtoTracks.empty()) {
    throw std::invalid_argument("Missing proto tracks input collections");
  }
  if (m_cfg.inputSourceLinks.empty()) {
    throw std::invalid_argument("Missing source links input collection");
  }
  if (m_cfg.outputTrackParameters.empty()) {
    throw std::invalid_argument("Missing track parameters output collection");
  }
  if (m_cfg.outputProtoTracks.empty()) {
    throw std::invalid_argument("Missing proto tracks output collections");
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

ActsExamples::ProcessCode ActsExamples::TrackParamsEstimationAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // Read input data
  const auto& protoTracks =
      ctx.eventStore.get<ProtoTrackContainer>(m_cfg.inputProtoTracks);
  // need source links to get the geometry identifer
  const auto& sourceLinks =
      ctx.eventStore.get<IndexSourceLinkContainer>(m_cfg.inputSourceLinks);

  TrackParametersContainer trackParameters;
  ProtoTrackContainer tracks;
  trackParameters.reserve(protoTracks.size());
  tracks.reserve(protoTracks.size());

  for (std::size_t itrack = 0; itrack < protoTracks.size(); ++itrack) {
    // The list of hits and the initial start parameters
    const auto& protoTrack = protoTracks[itrack];
    // Space points on the proto track
    std::vector<SimSpacePoint> spacePointsOnTrack;
    spacePointsOnTrack.reserve(protoTrack.size());
    // Loop over the hit index on the proto track
    for (const auto& hitIndex : protoTrack) {
      // Loop over the sets of space point container to find the space point
      for (const auto& isp : m_cfg.inputSpacePoints) {
        const auto& spacePoints =
            ctx.eventStore.get<SimSpacePointContainer>(isp);

        auto it =
            std::find_if(spacePoints.begin(), spacePoints.end(),
                         [&](const SimSpacePoint& spacePoint) {
                           return (spacePoint.measurementIndex() == hitIndex);
                         });
        if (it != spacePoints.end()) {
          spacePointsOnTrack.push_back(*it);
          break;
        }
      }
    }

    // At least three space points are required
    if (spacePointsOnTrack.size() < 3) {
      continue;
    }

    // Sort the space points
    std::sort(spacePointsOnTrack.begin(), spacePointsOnTrack.end(),
              [](const SimSpacePoint& lhs, const SimSpacePoint& rhs) {
                return lhs.r() < rhs.r();
              });

    // Create a three-spacepoints seed
    Acts::Seed<SimSpacePoint> seed(spacePointsOnTrack[0], spacePointsOnTrack[1],
                                   spacePointsOnTrack[2],
                                   spacePointsOnTrack[1].z());

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
    Acts::Vector3 field =
        getField(Acts::Vector3(bottomSP->x(), bottomSP->y(), bottomSP->z()));
    // Estimate the track parameters from seed
    auto optParams = Acts::estimateTrackParamsFromSeed(
        ctx.geoContext, seed.sp().begin(), seed.sp().end(), *surface, field,
        m_cfg.bFieldMin);
    if (not optParams.has_value()) {
      ACTS_WARNING("Estimation of track parameters for proto track "
                   << itrack << " failed.");
      continue;
    } else {
      const auto& params = optParams.value();
      double charge = std::copysign(1, params[Acts::eBoundQOverP]);
      trackParameters.emplace_back(surface->getSharedPtr(), params, charge,
                                   m_covariance);
      // Create a new proto track with only the three space points in the seed
      ProtoTrack track(3);
      for (size_t i = 0; i < 3; ++i) {
        track[i] = spacePointsOnTrack[i].measurementIndex();
      }
      tracks.emplace_back(track);
    }
  }

  ctx.eventStore.add(m_cfg.outputTrackParameters, std::move(trackParameters));
  ctx.eventStore.add(m_cfg.outputProtoTracks, std::move(tracks));
  return ActsExamples::ProcessCode::SUCCESS;
}
