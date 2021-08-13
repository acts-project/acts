// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"

#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <stdexcept>

ActsExamples::TrackFindingAlgorithm::TrackFindingAlgorithm(
    Config cfg, Acts::Logging::Level level)
    : ActsExamples::BareAlgorithm("TrackFindingAlgorithm", level),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurements input collection");
  }
  if (m_cfg.inputSourceLinks.empty()) {
    throw std::invalid_argument("Missing source links input collection");
  }
  if (m_cfg.inputInitialTrackParameters.empty()) {
    throw std::invalid_argument(
        "Missing initial track parameters input collection");
  }
  if (m_cfg.outputTrajectories.empty()) {
    throw std::invalid_argument("Missing trajectories output collection");
  }
}

ActsExamples::ProcessCode ActsExamples::TrackFindingAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // Read input data
  const auto& measurements =
      ctx.eventStore.get<MeasurementContainer>(m_cfg.inputMeasurements);
  const auto& sourceLinks =
      ctx.eventStore.get<IndexSourceLinkContainer>(m_cfg.inputSourceLinks);
  const auto& initialParameters = ctx.eventStore.get<TrackParametersContainer>(
      m_cfg.inputInitialTrackParameters);

  // Prepare the output data with MultiTrajectory
  TrajectoriesContainer trajectories;
  trajectories.reserve(initialParameters.size());

  // Construct a perigee surface as the target surface
  auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
      Acts::Vector3{0., 0., 0.});

  Acts::PropagatorPlainOptions pOptions;
  pOptions.maxSteps = 10000;

  // Set the CombinatorialKalmanFilter options
  ActsExamples::TrackFindingAlgorithm::TrackFinderOptions options(
      ctx.geoContext, ctx.magFieldContext, ctx.calibContext,
      IndexSourceLinkAccessor(), MeasurementCalibrator(measurements),
      Acts::MeasurementSelector(m_cfg.measurementSelectorCfg),
      Acts::LoggerWrapper{logger()}, pOptions, &(*pSurface));

  // Perform the track finding for all initial parameters
  ACTS_DEBUG("Invoke track finding with " << initialParameters.size()
                                          << " seeds.");
  auto results = m_cfg.findTracks(sourceLinks, initialParameters, options);

  // Compute shared hits from all the reconstructed tracks
  // Compute nSharedhits and Update ckf results
  // hit index -> list of multi traj indexes [traj, meas]
  // @todo put shared hit computation in ad-hoc method
  if (m_cfg.computeSharedHits) {
    std::vector<int> firstTrackOnTheHit(sourceLinks.size(), -1);
    std::vector<int> firstStateOnTheHit(sourceLinks.size(), -1);

    for (unsigned int iresult(0); iresult < results.size(); iresult++) {
      if (not results.at(iresult).ok())
        continue;

      auto& ckfResult = results.at(iresult).value();
      auto& measIndexes = ckfResult.lastMeasurementIndices;

      for (auto measIndex : measIndexes) {
        ckfResult.fittedStates.visitBackwards(
            measIndex, [&](const auto& state) {
              if (not state.typeFlags().test(
                      Acts::TrackStateFlag::MeasurementFlag))
                return;

              std::size_t hitIndex = state.uncalibrated().index();

              // Check if hit not already used
              if (firstTrackOnTheHit.at(hitIndex) == -1) {
                firstTrackOnTheHit.at(hitIndex) = iresult;
                firstStateOnTheHit.at(hitIndex) = state.index();
                return;
              }

              // if already used, control if first track state has been marked
              // as shared
              int indexFirstTrack = firstTrackOnTheHit.at(hitIndex);
              int indexFirstState = firstStateOnTheHit.at(hitIndex);
              if (not results.at(indexFirstTrack)
                          .value()
                          .fittedStates.getTrackState(indexFirstState)
                          .typeFlags()
                          .test(Acts::TrackStateFlag::SharedHitFlag))
                results.at(indexFirstTrack)
                    .value()
                    .fittedStates.getTrackState(indexFirstState)
                    .typeFlags()
                    .set(Acts::TrackStateFlag::SharedHitFlag);

              // Decorate this track
              results.at(iresult)
                  .value()
                  .fittedStates.getTrackState(state.index())
                  .typeFlags()
                  .set(Acts::TrackStateFlag::SharedHitFlag);
            });
      }
    }
  }

  // Loop over the track finding results for all initial parameters
  for (std::size_t iseed = 0; iseed < initialParameters.size(); ++iseed) {
    // The result for this seed
    auto& result = results[iseed];
    if (result.ok()) {
      // Get the track finding output object
      const auto& trackFindingOutput = result.value();
      // Create a Trajectories result struct
      trajectories.emplace_back(
          std::move(trackFindingOutput.fittedStates),
          std::move(trackFindingOutput.lastMeasurementIndices),
          std::move(trackFindingOutput.fittedParameters));
    } else {
      ACTS_WARNING("Track finding failed for seed " << iseed << " with error"
                                                    << result.error());
      // Track finding failed. Add an empty result so the output container has
      // the same number of entries as the input.
      trajectories.push_back(Trajectories());
    }
  }

  ACTS_DEBUG("Finalized track finding with " << trajectories.size()
                                             << " track candidates.");

  ctx.eventStore.add(m_cfg.outputTrajectories, std::move(trajectories));
  return ActsExamples::ProcessCode::SUCCESS;
}
