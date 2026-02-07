// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Alignment/AlignmentAlgorithm.hpp"

#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "ActsAlignment/Kernel/AlignmentMask.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/MeasurementCalibration.hpp"
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
  if (m_cfg.trackingGeometry == nullptr) {
    throw std::invalid_argument("Missing tracking geometry");
  }

  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputProtoTracks.initialize(m_cfg.inputProtoTracks);
  m_inputInitialTrackParameters.initialize(m_cfg.inputInitialTrackParameters);
  m_outputAlignmentParameters.initialize(m_cfg.outputAlignmentParameters);
}

ActsExamples::ProcessCode ActsExamples::AlignmentAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // Read input data
  const auto& measurements = m_inputMeasurements(ctx);
  const auto& protoTracks = m_inputProtoTracks(ctx);
  const auto& initialParameters = m_inputInitialTrackParameters(ctx);

  // Save contexts from the first event
  if (m_collectedSourceLinks.empty()) {
    m_savedGeoContext = ctx.geoContext;
    m_savedMagFieldContext = ctx.magFieldContext;
    m_savedCalibContext = ctx.calibContext;
    // Initialize measurements container (will be populated with all events)
    m_collectedMeasurements = std::make_shared<MeasurementContainer>();
  }

  // Consistency cross checks
  if (protoTracks.size() != initialParameters.size()) {
    ACTS_FATAL("Inconsistent number of proto tracks and parameters "
               << protoTracks.size() << " vs " << initialParameters.size());
    return ProcessCode::ABORT;
  }

  std::size_t numTracksUsed = protoTracks.size();
  if (m_cfg.maxNumTracks > 0 &&
      m_collectedSourceLinks.size() + numTracksUsed >
          static_cast<std::size_t>(m_cfg.maxNumTracks)) {
    numTracksUsed = static_cast<std::size_t>(m_cfg.maxNumTracks) -
                    m_collectedSourceLinks.size();
  }

  ACTS_DEBUG("Collecting " << numTracksUsed << " tracks from event "
                           << ctx.eventNumber);

  // Get the current size of collected measurements (for index offset)
  std::size_t measurementOffset = m_collectedMeasurements->size();
  ACTS_DEBUG("Current collected measurements size: "
             << measurementOffset
             << ", event measurements size: " << measurements.size());

  // Merge measurements from this event into the collected container
  // Create a mapping from old index to new index
  std::vector<std::size_t> indexMap(measurements.size());
  for (std::size_t oldIdx = 0; oldIdx < measurements.size(); ++oldIdx) {
    const ConstVariableBoundMeasurementProxy oldMeasurement =
        measurements.getMeasurement(oldIdx);

    // Copy measurement to collected container using copyMeasurement
    auto newMeasurementProxy =
        m_collectedMeasurements->copyMeasurement(oldMeasurement);
    indexMap[oldIdx] = newMeasurementProxy.index();
  }

  ACTS_DEBUG("Merged " << measurements.size() << " measurements. "
                       << "New collected measurements size: "
                       << m_collectedMeasurements->size());

  // Collect track data with remapped indices
  std::vector<IndexSourceLink> trackSourceLinks;
  for (std::size_t itrack = 0; itrack < numTracksUsed; ++itrack) {
    // The list of hits and the initial start parameters
    const auto& protoTrack = protoTracks[itrack];

    // Clear & reserve the right size
    trackSourceLinks.clear();
    trackSourceLinks.reserve(protoTrack.size());

    // Fill the source links with remapped indices
    for (auto measIndex : protoTrack) {
      const ConstVariableBoundMeasurementProxy measurement =
          measurements.getMeasurement(measIndex);
      // Use the remapped index instead of the original index
      std::size_t remappedIndex = indexMap[measIndex];
      IndexSourceLink sourceLink(measurement.geometryId(), remappedIndex);
      trackSourceLinks.push_back(sourceLink);
    }

    // Store this track's data
    m_collectedSourceLinks.push_back(trackSourceLinks);
    m_collectedInitialParameters.push_back(initialParameters[itrack]);
  }

  ACTS_INFO("Total collected tracks so far: " << m_collectedSourceLinks.size());

  return ActsExamples::ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::AlignmentAlgorithm::finalize() {
  ACTS_INFO("=============================================================");
  ACTS_INFO("Finalizing alignment with "
            << m_collectedSourceLinks.size()
            << " collected tracks from all events");
  ACTS_INFO("=============================================================");

  if (m_collectedSourceLinks.empty()) {
    ACTS_WARNING("No tracks collected for alignment!");
    return ProcessCode::SUCCESS;
  }

  // Prepare the output for alignment parameters
  AlignmentParameters alignedParameters;

  // Construct a perigee surface as the target surface for the fitter
  auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
      Acts::Vector3{0., 0., 0.});

  Acts::KalmanFitterExtensions<Acts::VectorMultiTrajectory> extensions;
  PassThroughCalibrator pcalibrator;
  MeasurementCalibratorAdapter calibrator(pcalibrator,
                                          *m_collectedMeasurements);
  extensions.calibrator.connect<&MeasurementCalibratorAdapter::calibrate>(
      &calibrator);
  Acts::GainMatrixUpdater kfUpdater;
  Acts::GainMatrixSmoother kfSmoother;
  extensions.updater.connect<
      &Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(
      &kfUpdater);
  extensions.smoother.connect<
      &Acts::GainMatrixSmoother::operator()<Acts::VectorMultiTrajectory>>(
      &kfSmoother);

  // Set up the surface accessor to get surfaces from IndexSourceLink
  IndexSourceLink::SurfaceAccessor surfaceAccessor{*m_cfg.trackingGeometry};
  extensions.surfaceAccessor
      .connect<&IndexSourceLink::SurfaceAccessor::operator()>(&surfaceAccessor);

  // Set the KalmanFitter options using saved contexts
  TrackFitterOptions kfOptions(
      m_savedGeoContext, m_savedMagFieldContext, m_savedCalibContext,
      extensions,
      Acts::PropagatorPlainOptions(m_savedGeoContext, m_savedMagFieldContext),
      &(*pSurface));

  // 转换 iterationState: std::map<uint, bitset<6>> -> std::map<uint,
  // AlignmentMask>
  std::map<unsigned int, ActsAlignment::AlignmentMask> alignIterationState;
  for (const auto& [iter, mask] : m_cfg.iterationState) {
    alignIterationState[iter] =
        static_cast<ActsAlignment::AlignmentMask>(mask.to_ulong());
  }

  // Set the alignment options (包括 iterationState)
  ActsAlignment::AlignmentOptions<TrackFitterOptions> alignOptions(
      kfOptions, m_cfg.alignedTransformUpdater, m_cfg.alignedDetElements,
      m_cfg.chi2ONdfCutOff, m_cfg.deltaChi2ONdfCutOff, m_cfg.maxNumIterations,
      alignIterationState);

  ACTS_INFO("Starting alignment iterations...");
  ACTS_INFO("  Number of tracks: " << m_collectedSourceLinks.size());
  ACTS_INFO(
      "  Total collected measurements: " << m_collectedMeasurements->size());
  ACTS_INFO("  Aligned detector elements: " << m_cfg.alignedDetElements.size());
  ACTS_INFO("  Max iterations: " << m_cfg.maxNumIterations);

  auto result = (*m_cfg.align)(m_collectedSourceLinks,
                               m_collectedInitialParameters, alignOptions);

  if (result.ok()) {
    const auto& alignOutput = result.value();
    alignedParameters = alignOutput.alignedParameters;
    ACTS_INFO("=============================================================");
    ACTS_INFO("Alignment finished successfully!");
    ACTS_INFO("  Final deltaChi2: " << result.value().deltaChi2);
    ACTS_INFO("  Aligned parameters: " << alignedParameters.size());
    ACTS_INFO("=============================================================");
  } else {
    ACTS_ERROR("Alignment failed with error: " << result.error());
    return ProcessCode::ABORT;
  }

  // TODO: Save alignment parameters to file
  // For now, they are stored in the alignedDetElements via
  // alignedTransformUpdater

  return ActsExamples::ProcessCode::SUCCESS;
}
