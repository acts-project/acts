// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFitting/RefittingAlgorithm.hpp"

#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <stdexcept>

ActsExamples::RefittingAlgorithm::RefittingAlgorithm(Config config,
                                                     Acts::Logging::Level level)
    : ActsExamples::IAlgorithm("TrackFittingAlgorithm", level),
      m_cfg(std::move(config)) {
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Missing input tracks collection");
  }
  if (m_cfg.outputTracks.empty()) {
    throw std::invalid_argument("Missing output tracks collection");
  }
}

ActsExamples::ProcessCode ActsExamples::RefittingAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // Read input data
  const auto& inputTracks =
      ctx.eventStore.get<ConstTrackContainer>(m_cfg.inputTracks);

  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
  TrackContainer tracks(trackContainer, trackStateContainer);

  // Perform the fit for each input track
  std::vector<Acts::SourceLink> trackSourceLinks;
  std::vector<const Acts::Surface*> surfSequence;
  AlreadyCalibratedCalibrator calibrator;

  auto itrack = 0ul;
  for (const auto track : inputTracks) {
    // Check if you are not in picking mode
    if (m_cfg.pickTrack > -1 and
        m_cfg.pickTrack != static_cast<int>(itrack++)) {
      continue;
    }

    TrackFitterFunction::GeneralFitterOptions options{
        ctx.geoContext, ctx.magFieldContext, ctx.calibContext,
        &track.referenceSurface(), Acts::PropagatorPlainOptions()};

    const Acts::BoundTrackParameters initialParams(
        track.referenceSurface().getSharedPtr(), track.parameters(),
        track.covariance());

    trackSourceLinks.clear();
    surfSequence.clear();
    calibrator.callibratedStates.clear();

    for (auto state : track.trackStates()) {
      surfSequence.push_back(&state.referenceSurface());

      if (not state.hasUncalibratedSourceLink()) {
        continue;
      }

      trackSourceLinks.push_back(state.uncalibratedSourceLink());
      const auto slIndex =
          state.uncalibratedSourceLink().get<IndexSourceLink>().index();
      auto [it, success] =
          calibrator.callibratedStates.insert({slIndex, std::move(state)});
      assert(success);
    }

    if (surfSequence.empty()) {
      ACTS_WARNING("Empty track " << itrack << " found.");
      continue;
    }

    ACTS_VERBOSE("Initial parameters: "
                 << initialParams.fourPosition(ctx.geoContext).transpose()
                 << " -> " << initialParams.unitDirection().transpose());

    ACTS_DEBUG("Invoke direct fitter for track " << itrack);
    auto result = (*m_cfg.fit)(trackSourceLinks, initialParams, options,
                               calibrator, surfSequence, tracks);

    if (result.ok()) {
      // Get the fit output object
      const auto& refittedTrack = result.value();
      if (refittedTrack.hasReferenceSurface()) {
        ACTS_VERBOSE("Refitted parameters for track " << itrack);
        ACTS_VERBOSE("  " << track.parameters().transpose());
      } else {
        ACTS_DEBUG("No refitted parameters for track " << itrack);
      }
    } else {
      ACTS_WARNING("Fit failed for track "
                   << itrack << " with error: " << result.error() << ", "
                   << result.error().message());
    }
  }

  std::stringstream ss;
  trackStateContainer->statistics().toStream(ss);
  ACTS_DEBUG(ss.str());

  ConstTrackContainer constTracks{
      std::make_shared<Acts::ConstVectorTrackContainer>(
          std::move(*trackContainer)),
      std::make_shared<Acts::ConstVectorMultiTrajectory>(
          std::move(*trackStateContainer))};

  ctx.eventStore.add(m_cfg.outputTracks, std::move(constTracks));
  return ActsExamples::ProcessCode::SUCCESS;
}
