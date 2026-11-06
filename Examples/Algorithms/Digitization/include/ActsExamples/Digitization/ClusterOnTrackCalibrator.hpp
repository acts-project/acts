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

namespace ActsExamples {

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
