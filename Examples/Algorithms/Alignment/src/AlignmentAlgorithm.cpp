// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Alignment/AlignmentAlgorithm.hpp"

#include "Acts/Definitions/Alignment.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/Utilities/detail/OstreamStateGuard.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"
#include "ActsAlignment/Kernel/AlignmentMask.hpp"
#include "ActsExamples/DetectorCommons/AlignmentContext.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/MeasurementCalibration.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <utility>
#include <vector>

namespace ActsExamples {
namespace {

constexpr std::array<const char*, Acts::eAlignmentSize> kAlignmentDofNames = {
    "dx", "dy", "dz", "rx", "ry", "rz"};

constexpr std::size_t kNotCopied = std::numeric_limits<std::size_t>::max();

// Copy the hits used by the first nTracks tracks, once each, and append the
// matching parameters. Returns the first measurement index outside source.
template <typename Tracks, typename Parameters, typename IndexOf>
std::optional<std::size_t> copyReferencedMeasurements(
    MeasurementContainer& dest,
    std::vector<std::vector<IndexSourceLink>>& destTracks,
    TrackParametersContainer& destParameters,
    const MeasurementContainer& source, const Tracks& tracks,
    Parameters&& parameters, std::size_t nTracks, IndexOf indexOf) {
  std::vector<std::size_t> copied(source.size(), kNotCopied);
  destTracks.reserve(destTracks.size() + nTracks);
  destParameters.reserve(destParameters.size() + nTracks);
  for (std::size_t itrack = 0; itrack < nTracks; ++itrack) {
    std::vector<IndexSourceLink> track;
    track.reserve(tracks[itrack].size());
    for (const auto& hit : tracks[itrack]) {
      const auto oldIndex = static_cast<std::size_t>(indexOf(hit));
      if (oldIndex >= source.size()) {
        return oldIndex;
      }
      const auto measurement = source.getMeasurement(oldIndex);
      std::size_t& index = copied[oldIndex];
      if (index == kNotCopied) {
        index = dest.copyMeasurement(measurement).index();
      }
      track.emplace_back(measurement.geometryId(), index);
    }
    destTracks.push_back(std::move(track));
    destParameters.push_back(std::move(parameters[itrack]));
  }
  return std::nullopt;
}

bool writeAlignmentIndex(const ActsAlignment::AlignmentResult& alignOutput,
                         std::ostream& out) {
  // protect stream from leaking modifiers
  Acts::detail::OstreamStateGuard guard(out);
  out << "# label  surface_index  dof  dof_name  geo_id\n";
  for (const auto& [surface, surfaceIndex] :
       ActsAlignment::orderedAlignmentSurfaces(alignOutput)) {
    if (surface == nullptr) {
      return false;
    }
    for (std::size_t dof = 0; dof < Acts::eAlignmentSize; ++dof) {
      const std::size_t label = surfaceIndex * Acts::eAlignmentSize + dof + 1;
      out << std::left << std::setw(8) << label << ' ' << std::setw(14)
          << surfaceIndex << ' ' << std::setw(4) << dof << ' ' << std::setw(8)
          << kAlignmentDofNames[dof] << ' ' << surface->geometryId() << '\n';
    }
  }
  return static_cast<bool>(out);
}

}  // namespace

AlignmentAlgorithm::AlignmentAlgorithm(
    Config cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("AlignmentAlgorithm", std::move(logger)),
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
  if (m_cfg.align == nullptr) {
    throw std::invalid_argument("Missing alignment function");
  }
  if (m_cfg.trackingGeometry == nullptr) {
    throw std::invalid_argument("Missing tracking geometry");
  }

  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputProtoTracks.initialize(m_cfg.inputProtoTracks);
  m_inputInitialTrackParameters.initialize(m_cfg.inputInitialTrackParameters);
}

ProcessCode AlignmentAlgorithm::execute(const AlgorithmContext& ctx) const {
  const auto& measurements = m_inputMeasurements(ctx);
  const auto& protoTracks = m_inputProtoTracks(ctx);
  const auto& initialParameters = m_inputInitialTrackParameters(ctx);

  if (protoTracks.size() != initialParameters.size()) {
    ACTS_FATAL("Inconsistent number of proto tracks and parameters "
               << protoTracks.size() << " vs " << initialParameters.size());
    return ProcessCode::ABORT;
  }

  // Each event keeps at most the cap. takeCollectedSample() then keeps the
  // event-number prefix.
  const std::size_t nTracks =
      m_cfg.maxNumTracks > 0
          ? std::min(protoTracks.size(),
                     static_cast<std::size_t>(m_cfg.maxNumTracks))
          : protoTracks.size();

  std::optional<std::size_t> badIndex;
  std::optional<CollectedChunk> chunk;
  if (nTracks > 0) {
    chunk.emplace();
    chunk->geoContext = ctx.recoGeoContext;
    chunk->magFieldContext = ctx.magFieldContext;
    chunk->calibContext = ctx.calibContext;
    badIndex = copyReferencedMeasurements(
        chunk->measurements, chunk->sourceLinks, chunk->parameters,
        measurements, protoTracks, initialParameters, nTracks,
        [](auto index) { return index; });
  }

  std::scoped_lock lock(m_collectionMutex);
  const IAlignmentStore* store = nullptr;
  if (const auto* alignment = ctx.recoGeoContext.maybeGet<AlignmentContext>()) {
    store = alignment->store;
  }
  if (m_alignmentStoreSeen && m_alignmentStore != store) {
    m_alignmentStoreMismatch = true;
    ACTS_ERROR("Event " << ctx.eventNumber
                        << " uses a different alignment store");
    return ProcessCode::ABORT;
  }
  m_alignmentStoreSeen = true;
  m_alignmentStore = store;
  if (badIndex.has_value()) {
    ACTS_ERROR("Measurement index "
               << *badIndex << " is outside the measurement container (size "
               << measurements.size() << ") in event " << ctx.eventNumber);
    return ProcessCode::ABORT;
  }
  if (!chunk.has_value()) {
    return ProcessCode::SUCCESS;
  }
  if (!m_collectedChunks.emplace(ctx.eventNumber, std::move(*chunk)).second) {
    ACTS_ERROR("Duplicate event number " << ctx.eventNumber);
    return ProcessCode::ABORT;
  }
  ACTS_DEBUG("Collecting " << nTracks << " tracks from event "
                           << ctx.eventNumber);
  return ProcessCode::SUCCESS;
}

AlignmentAlgorithm::CollectedChunk AlignmentAlgorithm::takeCollectedSample() {
  CollectedChunk sample;
  if (m_collectedChunks.empty()) {
    return sample;
  }
  const auto& first = m_collectedChunks.begin()->second;
  sample.geoContext = first.geoContext;
  sample.magFieldContext = first.magFieldContext;
  sample.calibContext = first.calibContext;

  const bool cap = m_cfg.maxNumTracks > 0;
  const auto limit = static_cast<std::size_t>(std::max(m_cfg.maxNumTracks, 0));
  std::size_t kept = 0;
  for (auto& [_, chunk] : m_collectedChunks) {
    if (cap && kept >= limit) {
      break;
    }
    const std::size_t nTracks =
        cap ? std::min(chunk.sourceLinks.size(), limit - kept)
            : chunk.sourceLinks.size();
    const auto badIndex = copyReferencedMeasurements(
        sample.measurements, sample.sourceLinks, sample.parameters,
        chunk.measurements, chunk.sourceLinks, chunk.parameters, nTracks,
        [](const IndexSourceLink& link) { return link.index(); });
    assert(!badIndex.has_value());
    (void)badIndex;
    kept += nTracks;
  }
  m_collectedChunks.clear();
  return sample;
}

ProcessCode AlignmentAlgorithm::initialize() {
  std::scoped_lock lock(m_collectionMutex);
  m_collectedChunks.clear();
  m_alignmentStore = nullptr;
  m_alignmentStoreSeen = false;
  m_alignmentStoreMismatch = false;
  return ProcessCode::SUCCESS;
}

ProcessCode AlignmentAlgorithm::finalize() {
  if (m_alignmentStoreMismatch) {
    ACTS_ERROR("Alignment store changed between events");
    return ProcessCode::ABORT;
  }

  const CollectedChunk sample = takeCollectedSample();
  if (sample.sourceLinks.empty()) {
    ACTS_WARNING("No tracks collected for alignment");
    return ProcessCode::SUCCESS;
  }

  auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
      Acts::Vector3{0., 0., 0.});

  Acts::KalmanFitterExtensions<Acts::VectorMultiTrajectory> extensions;
  PassThroughCalibrator pcalibrator;
  MeasurementCalibratorAdapter calibrator(pcalibrator, sample.measurements);
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

  IndexSourceLink::SurfaceAccessor surfaceAccessor{*m_cfg.trackingGeometry};
  extensions.surfaceAccessor
      .connect<&IndexSourceLink::SurfaceAccessor::operator()>(&surfaceAccessor);

  TrackFitterOptions kfOptions(
      sample.geoContext, sample.magFieldContext, sample.calibContext,
      extensions,
      Acts::PropagatorPlainOptions(sample.geoContext, sample.magFieldContext),
      &(*pSurface));

  std::map<unsigned int, ActsAlignment::AlignmentMask> iterationState;
  for (const auto& [iteration, mask] : m_cfg.iterationState) {
    iterationState.emplace(
        iteration, static_cast<ActsAlignment::AlignmentMask>(mask.to_ulong()));
  }
  ActsAlignment::AlignmentOptions<TrackFitterOptions> alignOptions(
      kfOptions, m_cfg.alignedTransformUpdater, m_cfg.alignedDetElements,
      m_cfg.chi2ONdfCutOff, m_cfg.deltaChi2ONdfCutOff, m_cfg.maxNumIterations,
      iterationState);

  ACTS_INFO("Alignment tracks: "
            << sample.sourceLinks.size()
            << ", measurements: " << sample.measurements.size()
            << ", detector elements: " << m_cfg.alignedDetElements.size());

  auto result =
      (*m_cfg.align)(sample.sourceLinks, sample.parameters, alignOptions);
  if (!result.ok()) {
    ACTS_ERROR("Alignment failed with " << result.error());
    return ProcessCode::ABORT;
  }

  const auto& alignOutput = result.value();
  ACTS_INFO("Alignment finished with deltaChi2 = "
            << alignOutput.deltaChi2
            << ", parameters: " << alignOutput.alignedParameters.size());
  if (!alignOutput.result.ok()) {
    ACTS_WARNING("Alignment is not converged ("
                 << alignOutput.result.error().message() << ")");
  }

  if (m_cfg.outputAlignmentFile.has_value()) {
    std::ofstream out(*m_cfg.outputAlignmentFile);
    if (!out) {
      ACTS_ERROR("Failed to open alignment result file: "
                 << *m_cfg.outputAlignmentFile);
      return ProcessCode::ABORT;
    }
    ActsAlignment::dumpResultAsMP2(alignOutput, out);
    if (!out) {
      ACTS_ERROR("Failed to write alignment result file: "
                 << *m_cfg.outputAlignmentFile);
      return ProcessCode::ABORT;
    }
    ACTS_INFO("Wrote alignment result file: " << *m_cfg.outputAlignmentFile);
  }
  if (m_cfg.outputAlignmentIndexFile.has_value()) {
    std::ofstream out(*m_cfg.outputAlignmentIndexFile);
    if (!out || !writeAlignmentIndex(alignOutput, out)) {
      ACTS_ERROR("Failed to write alignment index file: "
                 << *m_cfg.outputAlignmentIndexFile);
      return ProcessCode::ABORT;
    }
    ACTS_INFO(
        "Wrote alignment index file: " << *m_cfg.outputAlignmentIndexFile);
  }
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
