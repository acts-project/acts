// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/TrackExtrapolationAlgorithm.hpp"

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

  const Options options(ctx.geoContext, ctx.magFieldContext);

  auto trackBackend = std::make_shared<Acts::VectorTrackContainer>();
  auto stateBackend = std::make_shared<Acts::VectorMultiTrajectory>();
  TrackContainer extrapolated{trackBackend, stateBackend};
  extrapolated.ensureDynamicColumns(inputTracks);

  std::size_t nFailed = 0;

  for (const auto& track : inputTracks) {
    auto destination = extrapolated.makeTrack();
    destination.copyFromWithoutStates(track);

    // `TrackProxy::copyFrom` would copy the states with a hardcoded
    // `TrackStatePropMask::All`, which throws on the states of a seed track
    // that hold no parameters at all
    for (const auto& source : track.trackStatesReversed()) {
      auto state = destination.appendTrackState(source.getMask());
      state.copyFrom(source, source.getMask(), true);
    }
    destination.reverseTrackStates();

    const auto result = Acts::extrapolateTrackToReferenceSurface(
        destination, *m_cfg.targetSurface, propagator, options, m_cfg.strategy,
        logger());
    if (!result.ok()) {
      ACTS_DEBUG("Extrapolation of track " << track.index() << " failed with "
                                           << result.error());
      ++nFailed;
      extrapolated.removeTrack(destination.index());
    }
  }

  if (nFailed > 0) {
    ACTS_DEBUG("Dropped " << nFailed << " of " << inputTracks.size()
                          << " tracks that could not be extrapolated");
  }

  ConstTrackContainer outputTracks{
      std::make_shared<Acts::ConstVectorTrackContainer>(
          std::move(*trackBackend)),
      std::make_shared<Acts::ConstVectorMultiTrajectory>(
          std::move(*stateBackend))};

  ACTS_DEBUG("Extrapolated " << outputTracks.size() << " tracks");

  m_outputTracks(ctx, std::move(outputTracks));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
