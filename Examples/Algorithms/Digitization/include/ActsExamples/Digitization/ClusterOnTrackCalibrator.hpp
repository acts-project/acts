// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Digitization/DigitizationConfig.hpp"
#include "ActsExamples/EventData/MeasurementCalibration.hpp"

#include <cstddef>

namespace ActsExamples {

namespace detail {

/// The projection of a cluster onto one axis of the segmentation.
struct ClusterProjection {
  /// The edges of the first and the last bin of the cluster on this axis
  double firstLower = 0;
  double firstUpper = 0;
  double lastLower = 0;
  double lastUpper = 0;
  /// The summed activation of the cells in the first and the last bin
  double firstCharge = 0;
  double lastCharge = 0;
  /// The number of cells in the first and the last bin
  std::size_t firstCells = 0;
  std::size_t lastCells = 0;
  /// The widths of the bins before the first and after the last bin. Zero if
  /// the sensor ends there.
  double beforeWidth = 0;
  double afterWidth = 0;
};

/// The footprint of the predicted track in the sensor, along one axis.
struct ClusterFootprint {
  /// The extent of the path through the sensor on this axis
  double length = 0;
  /// The error of the extent, from the error of the predicted direction
  double lengthError = 0;
  /// The full length of the path through the sensor
  double pathLength = 0;
};

/// The calibrated position on one axis.
struct ClusterAxisCalibration {
  double position = 0;
  double variance = 0;
  /// False if the cluster cannot come from the predicted footprint
  bool consistent = false;
};

/// Calibrate the position of a cluster on one axis of the segmentation.
///
/// The track deposits a straight segment with a uniform charge density. The
/// digitization adds Gaussian noise to the charge of each cell that the
/// segment touches, and a cell fires if its charge is above the threshold.
/// The result is the mean and the variance of the track position, for a
/// uniform prior and this likelihood of the cluster:
///
/// - The first and the last bin fire. For analogue readout, their charge is
///   Gaussian around the charge of the predicted overlap instead.
/// - The touched bins before the first and after the last bin do not fire.
///
/// Without noise, the likelihood is 0 or 1, and the posterior is uniform in a
/// box. With noise, the posterior is integrated numerically. The error of the
/// predicted footprint length adds to the variance.
///
/// If no position gives the observed cluster, the cluster cannot come from
/// this footprint. The result is then uniform over the cluster.
///
/// @param projection The cluster on this axis
/// @param footprint The predicted track in the sensor on this axis
/// @param threshold The activation threshold of a cell
/// @param analogue Whether to use the charge of the cells
/// @param chargeNoise The Gaussian noise of the activation of one cell
ClusterAxisCalibration calibrateClusterAxis(const ClusterProjection& projection,
                                            const ClusterFootprint& footprint,
                                            double threshold, bool analogue,
                                            double chargeNoise);

}  // namespace detail

/// Calibrate geometric clusters with the direction of the predicted track.
///
/// The calibrator reads the sensor thickness, the segmentation, the threshold
/// and the readout mode from the geometric digitization configuration. It
/// supports plane surfaces with a Cartesian segmentation. All other
/// measurements, and measurements without track parameters, pass through.
///
/// Each axis is calibrated on its own. The cells of the cluster in one bin of
/// the axis count as one cell with the sum of their charge.
class ClusterOnTrackCalibrator final : public MeasurementCalibrator {
 public:
  struct Config {
    /// The digitization configuration that made the clusters
    DigiConfigContainer digitizationConfigs;
    /// The Gaussian noise of the activation of one cell. It must match the
    /// charge smearing of the digitization.
    double chargeNoise = 0;
    /// The minimum error of a calibrated position. Without charge noise and
    /// without a direction error, the error of a position can be zero.
    double minError = 1 * Acts::UnitConstants::um;
  };

  explicit ClusterOnTrackCalibrator(const Config& config);

  void calibrate(
      const MeasurementContainer& measurements,
      const ClusterContainer* clusters, const Acts::GeometryContext& gctx,
      const Acts::CalibrationContext& cctx, const Acts::SourceLink& sourceLink,
      Acts::VectorMultiTrajectory::TrackStateProxy& trackState) const override;

  bool needsClusters() const override { return true; }

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
