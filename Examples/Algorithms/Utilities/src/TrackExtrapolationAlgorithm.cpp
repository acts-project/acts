// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/TrackExtrapolationAlgorithm.hpp"

#include "Acts/Definitions/Direction.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/ActorList.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/SympyStepper.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Result.hpp"

#include <memory>
#include <stdexcept>
#include <utility>

namespace ActsExamples {

TrackExtrapolationAlgorithm::TrackExtrapolationAlgorithm(
    Config config, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("TrackExtrapolationAlgorithm", std::move(logger)),
      m_cfg(std::move(config)) {
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Missing input track collection");
  }
  if (m_cfg.outputTracks.empty()) {
    throw std::invalid_argument("Missing output track collection");
  }
  if (m_cfg.targetSurface == nullptr) {
    throw std::invalid_argument("Missing target surface");
  }
  if (m_cfg.trackingGeometry == nullptr) {
    throw std::invalid_argument("Missing tracking geometry");
  }
  if (m_cfg.magneticField == nullptr) {
    throw std::invalid_argument("Missing magnetic field");
  }

  m_inputTracks.initialize(m_cfg.inputTracks);
  m_outputTracks.initialize(m_cfg.outputTracks);
}

ProcessCode TrackExtrapolationAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  const ConstTrackContainer& inputTracks = m_inputTracks(ctx);

  using Propagator = Acts::Propagator<Acts::SympyStepper, Acts::Navigator>;
  using Options = Propagator::Options<
      Acts::ActorList<Acts::MaterialInteractor, Acts::EndOfWorldReached>>;

  const Propagator propagator(
      Acts::SympyStepper(m_cfg.magneticField),
      Acts::Navigator({m_cfg.trackingGeometry},
                      logger().cloneWithSuffix("Navigator")),
      logger().cloneWithSuffix("Propagator"));

  Options options(ctx.geoContext, ctx.magFieldContext);
  options.constrainToVolumeIds = m_cfg.constrainToVolumeIds;
  options.endOfWorldVolumeIds = m_cfg.endOfWorldVolumeIds;

  // `Acts::extrapolateTrackToReferenceSurface` split up, so the states are read
  // off the input track and only the parameters are written to the output one
  auto extrapolate = [&](const ConstTrackProxy& track)
      -> Acts::Result<Acts::BoundTrackParameters> {
    auto findResult = Acts::findTrackStateForExtrapolation(
        ctx.geoContext, track, *m_cfg.targetSurface, m_cfg.strategy, logger());
    if (!findResult.ok()) {
      return findResult.error();
    }
    const auto& [trackState, distance] = *findResult;

    Options trackOptions = options;
    trackOptions.direction =
        Acts::Direction::fromScalarZeroAsPositive(distance);

    auto propagateResult =
        propagator.propagate<Options, Acts::ForcedSurfaceReached>(
            track.createParametersFromState(trackState), *m_cfg.targetSurface,
            trackOptions);
    if (!propagateResult.ok()) {
      return propagateResult.error();
    }

    return propagateResult->endParameters.value();
  };

  // The tip and stem indices point into the input state backend, which the
  // output container takes over below. The empty backend here is only because
  // `Acts::TrackContainer` cannot pair a mutable track backend with a
  // read-only state one.
  auto trackBackend = std::make_shared<Acts::VectorTrackContainer>();
  TrackContainer extrapolated{trackBackend,
                              std::make_shared<Acts::VectorMultiTrajectory>()};
  extrapolated.ensureDynamicColumns(inputTracks);

  std::size_t nFailed = 0;

  for (const auto& track : inputTracks) {
    const auto result = extrapolate(track);
    if (!result.ok()) {
      ACTS_DEBUG("Extrapolation of track " << track.index() << " failed with "
                                           << result.error());
      ++nFailed;
      continue;
    }

    auto destination = extrapolated.makeTrack();
    destination.copyFromShallow(track);
    destination.setReferenceSurface(m_cfg.targetSurface);
    destination.parameters() = result->parameters();
    destination.covariance() = result->covariance().value();
  }

  m_nTotalTracks += inputTracks.size();
  m_nFailedTracks += nFailed;

  if (nFailed > 0) {
    ACTS_DEBUG("Dropped " << nFailed << " of " << inputTracks.size()
                          << " tracks that could not be extrapolated");
  }

  ConstTrackContainer outputTracks{
      std::make_shared<Acts::ConstVectorTrackContainer>(
          std::move(*trackBackend)),
      inputTracks.trackStateContainerHolder()};

  ACTS_DEBUG("Extrapolated " << outputTracks.size() << " tracks");

  m_outputTracks(ctx, std::move(outputTracks));

  return ProcessCode::SUCCESS;
}

ProcessCode TrackExtrapolationAlgorithm::finalize() {
  ACTS_INFO("TrackExtrapolationAlgorithm statistics:");
  ACTS_INFO("- total tracks: " << m_nTotalTracks);
  ACTS_INFO("- dropped tracks: " << m_nFailedTracks);
  if (m_nTotalTracks > 0) {
    ACTS_INFO("- drop ratio: " << static_cast<double>(m_nFailedTracks) /
                                      m_nTotalTracks);
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
