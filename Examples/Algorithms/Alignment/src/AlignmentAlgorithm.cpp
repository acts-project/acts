// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Alignment/AlignmentAlgorithm.hpp"

#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

ActsExamples::AlignmentAlgorithm::AlignmentAlgorithm(Config cfg,
                                                     Acts::Logging::Level lvl)
    : ActsExamples::IAlgorithm("AlignmentAlgorithm", lvl),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument("Missing input measurement collection");
  }
  if (m_cfg.inputSourceLinks.empty()) {
    throw std::invalid_argument("Missing input source links collection");
  }
  if (m_cfg.inputProtoTracks.empty()) {
    throw std::invalid_argument("Missing input proto tracks collection");
  }
  if (m_cfg.inputInitialTrackParameters.empty()) {
    throw std::invalid_argument(
        "Missing input initial track parameters collection");
  }
  if (m_cfg.outputAlignmentParameters.empty()) {
    throw std::invalid_argument(
        "Missing output alignment parameters collection");
  }

  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);
  m_inputProtoTracks.initialize(m_cfg.inputProtoTracks);
  m_inputInitialTrackParameters.initialize(m_cfg.inputInitialTrackParameters);
  m_outputAlignmentParameters.initialize(m_cfg.outputAlignmentParameters);
}

ActsExamples::ProcessCode ActsExamples::AlignmentAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // Read input data
  const auto& measurements = m_inputMeasurements(ctx);
  const auto& sourceLinks = m_inputSourceLinks(ctx);
  const auto& protoTracks = m_inputProtoTracks(ctx);
  const auto& initialParameters = m_inputInitialTrackParameters(ctx);

  // Consistency cross checks
  if (protoTracks.size() != initialParameters.size()) {
    ACTS_FATAL("Inconsistent number of proto tracks and parameters "
               << protoTracks.size() << " vs " << initialParameters.size());
    return ProcessCode::ABORT;
  }

  size_t numTracksUsed = protoTracks.size();
  if (m_cfg.maxNumTracks > 0 and
      m_cfg.maxNumTracks < static_cast<int>(protoTracks.size())) {
    numTracksUsed = m_cfg.maxNumTracks;
  }

  // Prepare the input track collection
  std::vector<std::vector<IndexSourceLink>> sourceLinkTrackContainer;
  sourceLinkTrackContainer.reserve(numTracksUsed);
  std::vector<IndexSourceLink> trackSourceLinks;
  for (std::size_t itrack = 0; itrack < numTracksUsed; ++itrack) {
    // The list of hits and the initial start parameters
    const auto& protoTrack = protoTracks[itrack];

    // Clear & reserve the right size
    trackSourceLinks.clear();
    trackSourceLinks.reserve(protoTrack.size());

    // Fill the source links via their indices from the container
    for (auto hitIndex : protoTrack) {
      auto sourceLink = sourceLinks.nth(hitIndex);
      if (sourceLink == sourceLinks.end()) {
        ACTS_FATAL("Proto track " << itrack << " contains invalid hit index"
                                  << hitIndex);
        return ProcessCode::ABORT;
      }
      trackSourceLinks.push_back(*sourceLink);
    }
    sourceLinkTrackContainer.push_back(trackSourceLinks);
  }

  // Prepare the output for alignment parameters
  AlignmentParameters alignedParameters;

  // Construct a perigee surface as the target surface for the fitter
  auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
      Acts::Vector3{0., 0., 0.});

  Acts::KalmanFitterExtensions<Acts::VectorMultiTrajectory> extensions;
  MeasurementCalibrator calibrator{measurements};
  extensions.calibrator.connect<&MeasurementCalibrator::calibrate>(&calibrator);
  Acts::GainMatrixUpdater kfUpdater;
  Acts::GainMatrixSmoother kfSmoother;
  extensions.updater.connect<
      &Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(
      &kfUpdater);
  extensions.smoother.connect<
      &Acts::GainMatrixSmoother::operator()<Acts::VectorMultiTrajectory>>(
      &kfSmoother);

  // Set the KalmanFitter options
  TrackFitterOptions kfOptions(ctx.geoContext, ctx.magFieldContext,
                               ctx.calibContext, extensions,
                               Acts::PropagatorPlainOptions(), &(*pSurface));

  // Set the alignment options
  ActsAlignment::AlignmentOptions<TrackFitterOptions> alignOptions(
      kfOptions, m_cfg.alignedTransformUpdater, m_cfg.alignedDetElements,
      m_cfg.chi2ONdfCutOff, m_cfg.deltaChi2ONdfCutOff, m_cfg.maxNumIterations);

  ACTS_DEBUG("Invoke track-based alignment with " << numTracksUsed
                                                  << " input tracks");
  auto result =
      (*m_cfg.align)(sourceLinkTrackContainer, initialParameters, alignOptions);
  if (result.ok()) {
    const auto& alignOutput = result.value();
    alignedParameters = alignOutput.alignedParameters;
    ACTS_VERBOSE(
        "Alignment finished with deltaChi2 = " << result.value().deltaChi2);
  } else {
    ACTS_WARNING("Alignment failed with " << result.error());
  }

  // add alignment parameters to event store
  m_outputAlignmentParameters(ctx, std::move(alignedParameters));
  return ActsExamples::ProcessCode::SUCCESS;
}
