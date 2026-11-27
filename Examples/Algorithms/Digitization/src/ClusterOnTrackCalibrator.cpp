// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/ClusterOnTrackCalibrator.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "ActsExamples/Digitization/ClusterAxisCalibration.hpp"
#include "ActsExamples/Digitization/Smearers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>

namespace ActsExamples {

namespace {

/// The noise of one cell from the charge smearing of the digitization, or
/// nothing if the smearing is not Gaussian
std::optional<double> chargeNoise(const GeometricConfig& geoCfg) {
  const auto& smearer = geoCfg.chargeSmearer;
  if (!smearer || smearer.target<Digitization::Exact>() != nullptr) {
    return 0.;
  }
  if (const auto* gauss = smearer.target<Digitization::Gauss>();
      gauss != nullptr) {
    return gauss->sigma;
  }
  return std::nullopt;
}

/// Project a cluster onto one axis of the segmentation
ClusterAxisProjection projectCluster(const Cluster& cluster,
                                     const Acts::IAxis& axis,
                                     std::size_t axisIndex) {
  std::size_t first = std::numeric_limits<std::size_t>::max();
  std::size_t last = 0;
  for (const auto& cell : cluster.channels) {
    first = std::min<std::size_t>(first, cell.bin[axisIndex]);
    last = std::max<std::size_t>(last, cell.bin[axisIndex]);
  }
  ClusterAxisProjection projection;
  // Cell bins are zero-based while the axis bin indices start at one
  projection.firstLower = axis.getBinLowerBound(first + 1);
  projection.firstUpper = axis.getBinUpperBound(first + 1);
  projection.lastLower = axis.getBinLowerBound(last + 1);
  projection.lastUpper = axis.getBinUpperBound(last + 1);
  if (first > 0) {
    projection.beforeWidth = axis.getBinWidth(first);
  }
  if (last + 1 < axis.getNBins()) {
    projection.afterWidth = axis.getBinWidth(last + 2);
  }
  projection.firstCells = 0;
  projection.lastCells = 0;
  for (const auto& cell : cluster.channels) {
    if (cell.bin[axisIndex] == first) {
      projection.firstCharge += cell.activation;
      ++projection.firstCells;
    }
    if (cell.bin[axisIndex] == last) {
      projection.lastCharge += cell.activation;
      ++projection.lastCells;
    }
  }
  return projection;
}

}  // namespace

ClusterOnTrackCalibrator::ClusterOnTrackCalibrator(const Config& config)
    : m_cfg(config) {
  if (!std::isfinite(m_cfg.minError) || m_cfg.minError < 0) {
    throw std::invalid_argument("Minimum error must be finite and nonnegative");
  }
  for (const auto& digiCfg : m_cfg.digitizationConfigs) {
    const GeometricConfig& geoCfg = digiCfg.geometricDigiConfig;
    if (geoCfg.indices.empty()) {
      continue;
    }
    const auto noise = chargeNoise(geoCfg);
    if (!noise.has_value() || !std::isfinite(*noise) || *noise < 0) {
      throw std::invalid_argument(
          "ClusterOnTrackCalibrator needs Exact or Gauss charge smearing");
    }
  }
}

void ClusterOnTrackCalibrator::calibrate(
    const MeasurementContainer& measurements, const ClusterContainer* clusters,
    const Acts::GeometryContext& gctx, const Acts::CalibrationContext& cctx,
    const Acts::SourceLink& sourceLink,
    Acts::VectorMultiTrajectory::TrackStateProxy& trackState) const {
  if (clusters == nullptr) {
    throw std::invalid_argument("ClusterOnTrackCalibrator needs clusters");
  }

  // Start from the measurement as it is
  PassThroughCalibrator{}.calibrate(measurements, clusters, gctx, cctx,
                                    sourceLink, trackState);

  // The fitters allocate the parameters that they calibrate with
  if (!trackState.hasPredicted() && !trackState.hasFiltered() &&
      !trackState.hasSmoothed()) {
    return;
  }
  const auto parameters = trackState.parameters();

  const Acts::Surface& surface = trackState.referenceSurface();
  if (surface.type() != Acts::Surface::Plane) {
    return;
  }
  auto digiIt = m_cfg.digitizationConfigs.find(surface.geometryId());
  if (digiIt == m_cfg.digitizationConfigs.end()) {
    return;
  }
  const GeometricConfig& geoCfg = digiIt->geometricDigiConfig;
  if (geoCfg.indices.empty() || geoCfg.segmentation == nullptr ||
      geoCfg.thickness <= 0) {
    return;
  }

  const std::size_t index = sourceLink.get<IndexSourceLink>().index();
  assert(index < clusters->size() && "Cluster index out of bounds");
  const Cluster& cluster = clusters->at(index);
  if (cluster.channels.empty()) {
    return;
  }

  // The predicted direction in the local frame. Its error is not used: the
  // fit can inflate it far above the actual error.
  const Acts::Vector3 dir =
      surface.localToGlobalTransform(gctx).linear().transpose() *
      Acts::makeDirectionFromPhiTheta(parameters[Acts::eBoundPhi],
                                      parameters[Acts::eBoundTheta]);
  // The clamp keeps the divisions below safe, because the compiler can
  // evaluate them before the check
  const double normal = std::max(std::abs(dir.z()), 1e-12);
  if (std::abs(dir.z()) < 1e-12) {
    return;
  }

  ClusterAxisReadout readout;
  readout.threshold = geoCfg.threshold;
  // The constructor checked the charge smearing
  readout.noise = chargeNoise(geoCfg).value_or(0.);
  readout.analogue = !geoCfg.digital;
  // The geometric digitization uses the path length as the charge
  readout.chargePerLength = 1;

  const auto measurement = measurements.getMeasurement(index);
  const std::size_t size = measurement.size();
  Acts::visit_measurement(size, [&](auto N) {
    constexpr std::size_t kSize = decltype(N)::value;
    Acts::Vector<kSize> calibrated = trackState.template calibrated<kSize>();
    Acts::SquareMatrix<kSize> calibratedCovariance =
        trackState.template calibratedCovariance<kSize>();

    for (Acts::BoundIndices bIndex : geoCfg.indices) {
      if ((bIndex != Acts::eBoundLoc0 && bIndex != Acts::eBoundLoc1) ||
          !measurement.contains(bIndex) ||
          static_cast<std::size_t>(bIndex) >= geoCfg.segmentation->getNAxes()) {
        continue;
      }
      const auto axisIndex = static_cast<std::size_t>(bIndex);
      const auto& axis = geoCfg.segmentation->getAxis(bIndex);
      const auto expected = bIndex == Acts::eBoundLoc0
                                ? Acts::AxisDirection::AxisX
                                : Acts::AxisDirection::AxisY;
      if (axis.getDirection() != expected) {
        continue;
      }

      const ClusterAxisProjection projection =
          projectCluster(cluster, axis, axisIndex);

      // The footprint of the predicted track on this axis
      ClusterAxisFootprint footprint;
      footprint.length = geoCfg.thickness * std::abs(dir[bIndex]) / normal;
      footprint.pathLength = geoCfg.thickness / normal;

      const ClusterAxisPosterior axisCalibration =
          calibrateClusterAxis(projection, footprint, readout);

      // The digitization gives a diagonal covariance, so this removes no
      // correlation
      const std::size_t i = measurement.indexOf(bIndex);
      calibrated[i] = axisCalibration.mean;
      calibratedCovariance.row(i).setZero();
      calibratedCovariance.col(i).setZero();
      calibratedCovariance(i, i) =
          std::max(axisCalibration.variance, m_cfg.minError * m_cfg.minError);
    }

    trackState.template calibrated<kSize>() = calibrated;
    trackState.template calibratedCovariance<kSize>() = calibratedCovariance;
  });
}

}  // namespace ActsExamples
