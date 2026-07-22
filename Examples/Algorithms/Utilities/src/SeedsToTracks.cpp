// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/SeedsToTracks.hpp"

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"

#include <sstream>
#include <stdexcept>
#include <utility>

namespace ActsExamples {

SeedsToTracks::SeedsToTracks(Config cfg,
                             std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("SeedsToTracks", std::move(logger)), m_cfg(std::move(cfg)) {
  m_inputSeeds.initialize(m_cfg.inputSeeds);
  m_inputTrackParameters.maybeInitialize(m_cfg.inputTrackParameters);
  m_outputTracks.initialize(m_cfg.outputTracks);
}

ProcessCode SeedsToTracks::execute(const AlgorithmContext& ctx) const {
  const SeedContainer& seeds = m_inputSeeds(ctx);
  ACTS_DEBUG("Received " << seeds.size() << " seeds");

  const TrackParametersContainer* trackParameters = nullptr;
  if (m_inputTrackParameters.isInitialized()) {
    trackParameters = &m_inputTrackParameters(ctx);

    if (trackParameters->size() != seeds.size()) {
      throw std::runtime_error(
          "Number of seeds and track parameters do not match");
    }
  }

  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto mtj = std::make_shared<Acts::VectorMultiTrajectory>();
  TrackContainer tracks(trackContainer, mtj);

  for (std::size_t i = 0; i < seeds.size(); ++i) {
    const auto seed = seeds.at(i);

    auto track = tracks.makeTrack();
    std::uint32_t nMeasurements = 0;

    for (const auto& sp : seed.spacePoints()) {
      for (const auto& sourceLink : sp.sourceLinks()) {
        // the estimated parameters are stored as predicted parameters on the
        // innermost track state, i.e. their reference surface
        const bool attachParameters =
            nMeasurements == 0 && trackParameters != nullptr;

        auto trackStateProxy = track.appendTrackState(
            attachParameters ? Acts::TrackStatePropMask::Predicted
                             : Acts::TrackStatePropMask::None);
        trackStateProxy.typeFlags().setIsMeasurement();
        trackStateProxy.setUncalibratedSourceLink(Acts::SourceLink(sourceLink));
        ++nMeasurements;

        if (m_cfg.trackingGeometry != nullptr) {
          const Acts::GeometryIdentifier geoId =
              sourceLink.get<IndexSourceLink>().geometryId();
          const Acts::Surface* surface =
              m_cfg.trackingGeometry->findSurface(geoId);
          if (surface == nullptr) {
            std::ostringstream oss;
            oss << "No surface found for source-link geometry id " << geoId;
            throw std::runtime_error(oss.str());
          }
          trackStateProxy.setReferenceSurface(surface->getSharedPtr());
        }

        if (attachParameters) {
          const auto& trackParams = trackParameters->at(i);
          trackStateProxy.predicted() = trackParams.parameters();
          trackStateProxy.predictedCovariance() =
              trackParams.covariance().value_or(Acts::BoundMatrix::Zero());
        }
      }
    }

    track.nMeasurements() = nMeasurements;
    track.nHoles() = 0;
    track.nOutliers() = 0;

    if (trackParameters != nullptr) {
      const auto& trackParams = trackParameters->at(i);

      track.setReferenceSurface(trackParams.referenceSurface().getSharedPtr());
      track.parameters() = trackParams.parameters();
      if (trackParams.covariance().has_value()) {
        track.covariance() = *trackParams.covariance();
      }
    }
  }

  ConstTrackContainer constTracks{
      std::make_shared<Acts::ConstVectorTrackContainer>(
          std::move(*trackContainer)),
      std::make_shared<Acts::ConstVectorMultiTrajectory>(std::move(*mtj))};

  ACTS_DEBUG("Produced " << constTracks.size() << " tracks");

  m_outputTracks(ctx, std::move(constTracks));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
