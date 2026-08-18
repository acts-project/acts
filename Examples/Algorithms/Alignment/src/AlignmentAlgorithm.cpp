// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Alignment/AlignmentAlgorithm.hpp"

#include "Acts/Definitions/Alignment.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "ActsAlignment/Kernel/AlignmentMask.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/MeasurementCalibration.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace ActsExamples {
namespace {

constexpr std::array<const char*, Acts::eAlignmentSize> kAlignmentDofNames = {
    "dx", "dy", "dz", "rx", "ry", "rz"};

std::string geoIdToString(Acts::GeometryIdentifier id) {
  std::ostringstream os;
  os << "vol=" << id.volume() << "|lay=" << id.layer()
     << "|sen=" << id.sensitive();
  if (id.extra() != 0u) {
    os << "|ext=" << id.extra();
  }
  return os.str();
}

/// Local [dx, dy, dz, rx, ry, rz] taking ``from`` to ``to``
/// (same convention as Alignment::updateAlignmentParameters / ODD script).
Acts::AlignmentVector localDeltasBetweenTransforms(const Acts::Transform3& from,
                                                   const Acts::Transform3& to) {
  const Acts::RotationMatrix3 R0 = from.rotation();
  const Acts::Vector3 t0 = from.translation();
  const Acts::RotationMatrix3 R1 = to.rotation();
  const Acts::Vector3 t1 = to.translation();

  Acts::AlignmentVector out = Acts::AlignmentVector::Zero();
  out.segment<3>(Acts::eAlignmentCenter0) = R0.transpose() * (t1 - t0);

  const Acts::RotationMatrix3 Rdelta = R0.transpose() * R1;
  const double ry = std::asin(std::clamp(-Rdelta(2, 0), -1.0, 1.0));
  const double cy = std::cos(ry);
  double rx = 0.;
  double rz = 0.;
  if (std::abs(cy) > 1e-8) {
    rx = std::atan2(Rdelta(2, 1), Rdelta(2, 2));
    rz = std::atan2(Rdelta(1, 0), Rdelta(0, 0));
  } else {
    rx = std::atan2(-Rdelta(1, 2), Rdelta(1, 1));
    rz = 0.;
  }
  out(Acts::eAlignmentRotation0) = rx;
  out(Acts::eAlignmentRotation1) = ry;
  out(Acts::eAlignmentRotation2) = rz;
  return out;
}

double covarianceError(const Acts::DynamicMatrix& cov, std::size_t row) {
  if (cov.rows() <= static_cast<Eigen::Index>(row) ||
      cov.cols() <= static_cast<Eigen::Index>(row)) {
    return 0.;
  }
  const double var =
      cov(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(row));
  if (!(var > 0.) || !std::isfinite(var)) {
    return 0.;
  }
  return std::sqrt(var);
}

}  // namespace

AlignmentAlgorithm::AlignmentAlgorithm(
    Config cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("AlignmentAlgorithm", std::move(logger)),
      m_cfg(std::move(cfg)),
      m_savedGeoContext(Acts::GeometryContext::dangerouslyDefaultConstruct()) {
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

ProcessCode AlignmentAlgorithm::execute(const AlgorithmContext& ctx) const {
  // Per-event whiteboard reads are thread-safe; shared collectors are not.
  const auto& measurements = m_inputMeasurements(ctx);
  const auto& protoTracks = m_inputProtoTracks(ctx);
  const auto& initialParameters = m_inputInitialTrackParameters(ctx);

  // Consistency cross checks (event-local)
  if (protoTracks.size() != initialParameters.size()) {
    ACTS_FATAL("Inconsistent number of proto tracks and parameters "
               << protoTracks.size() << " vs " << initialParameters.size());
    return ProcessCode::ABORT;
  }

  // Serialize merges into the shared cross-event collectors. Without this,
  // Sequencer(numThreads!=1) races on push_back / copyMeasurement (double
  // free).
  std::scoped_lock lock(m_collectionMutex);

  // Save contexts from the first event that reaches the collector
  if (m_collectedSourceLinks.empty()) {
    m_savedGeoContext = ctx.geoContext;
    m_savedMagFieldContext = ctx.magFieldContext;
    m_savedCalibContext = ctx.calibContext;
    m_collectedMeasurements = std::make_shared<MeasurementContainer>();
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

  return ProcessCode::SUCCESS;
}

void AlignmentAlgorithm::writeAlignmentTextFiles(
    const ActsAlignment::AlignmentResult& alignOutput,
    const std::vector<Acts::Transform3>& startTransforms) const {
  const std::size_t nSurfaces = m_cfg.alignedDetElements.size();
  if (nSurfaces == 0u) {
    ACTS_WARNING("No aligned detector elements; skip alignment text output");
    return;
  }
  if (startTransforms.size() != nSurfaces) {
    ACTS_ERROR("Start transform snapshot size mismatch");
    return;
  }

  const auto& cov = alignOutput.alignmentCovariance;

  if (!m_cfg.outputAlignmentFile.empty()) {
    std::ofstream out(m_cfg.outputAlignmentFile);
    if (!out) {
      ACTS_ERROR("Failed to open alignment result file: "
                 << m_cfg.outputAlignmentFile);
    } else {
      for (std::size_t sidx = 0; sidx < nSurfaces; ++sidx) {
        auto* det = m_cfg.alignedDetElements[sidx];
        Acts::Transform3 finalTrf = startTransforms[sidx];
        if (auto it = alignOutput.alignedParameters.find(det);
            it != alignOutput.alignedParameters.end()) {
          finalTrf = it->second;
        }
        const Acts::AlignmentVector deltas =
            localDeltasBetweenTransforms(startTransforms[sidx], finalTrf);

        for (std::size_t dof = 0; dof < Acts::eAlignmentSize; ++dof) {
          const std::size_t label = sidx * Acts::eAlignmentSize + dof + 1;
          const std::size_t row = sidx * Acts::eAlignmentSize + dof;
          const double value = deltas(static_cast<Eigen::Index>(dof));
          const double error = covarianceError(cov, row);
          out << std::setw(8) << label << "   " << std::setw(12) << value
              << "   " << std::setw(4) << 0.0 << "   " << std::setw(12) << value
              << "   " << std::setw(12) << error << "    99 \n";
        }
      }
      ACTS_INFO("Wrote alignment result text: " << m_cfg.outputAlignmentFile);
    }
  }

  if (!m_cfg.outputAlignmentIndexFile.empty()) {
    std::ofstream out(m_cfg.outputAlignmentIndexFile);
    if (!out) {
      ACTS_ERROR("Failed to open alignment index file: "
                 << m_cfg.outputAlignmentIndexFile);
    } else {
      out << "# label  surface_index  dof  dof_name  geo_id\n";
      for (std::size_t sidx = 0; sidx < nSurfaces; ++sidx) {
        const auto& surface = m_cfg.alignedDetElements[sidx]->surface();
        const std::string geoId = geoIdToString(surface.geometryId());
        for (std::size_t dof = 0; dof < Acts::eAlignmentSize; ++dof) {
          const std::size_t label = sidx * Acts::eAlignmentSize + dof + 1;
          out << std::left << std::setw(8) << label << ' ' << std::setw(14)
              << sidx << ' ' << std::setw(4) << dof << ' ' << std::setw(8)
              << kAlignmentDofNames[dof] << ' ' << geoId << '\n';
        }
      }
      ACTS_INFO(
          "Wrote alignment index map text: " << m_cfg.outputAlignmentIndexFile);
    }
  }
}

ProcessCode AlignmentAlgorithm::finalize() {
  ACTS_INFO("=============================================================");
  ACTS_INFO("Finalizing alignment with "
            << m_collectedSourceLinks.size()
            << " collected tracks from all events");
  ACTS_INFO("=============================================================");

  if (m_collectedSourceLinks.empty()) {
    ACTS_WARNING("No tracks collected for alignment!");
    return ProcessCode::SUCCESS;
  }

  // Snapshot poses before fit so result deltas are relative to the start
  // geometry (misaligned input), matching the ODD steering script convention.
  std::vector<Acts::Transform3> startTransforms;
  startTransforms.reserve(m_cfg.alignedDetElements.size());
  for (auto* det : m_cfg.alignedDetElements) {
    startTransforms.push_back(
        det->surface().localToGlobalTransform(m_savedGeoContext));
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

  // Convert iterationState: std::map<uint, bitset<6>> -> std::map<uint,
  // AlignmentMask>
  std::map<unsigned int, ActsAlignment::AlignmentMask> alignIterationState;
  for (const auto& [iter, mask] : m_cfg.iterationState) {
    alignIterationState[iter] =
        static_cast<ActsAlignment::AlignmentMask>(mask.to_ulong());
  }

  // Set the alignment options (including iterationState)
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

  if (!result.ok()) {
    ACTS_ERROR("Alignment failed with error: " << result.error());
    return ProcessCode::ABORT;
  }

  const auto& alignOutput = result.value();
  alignedParameters = alignOutput.alignedParameters;
  ACTS_INFO("=============================================================");
  ACTS_INFO("Alignment finished successfully!");
  ACTS_INFO("  Final deltaChi2: " << alignOutput.deltaChi2);
  ACTS_INFO("  Aligned parameters: " << alignedParameters.size());
  ACTS_INFO("  Alignment DoF: " << alignOutput.alignmentDof);
  ACTS_INFO("=============================================================");

  // Persist MP-style result (+ index map) from finalize (R8 / #5432).
  writeAlignmentTextFiles(alignOutput, startTransforms);

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
