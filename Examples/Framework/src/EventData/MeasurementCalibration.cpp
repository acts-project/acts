// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <ActsExamples/EventData/MeasurementCalibration.hpp>

void ActsExamples::PassThroughCalibrator::calibrate(
    const MeasurementContainer& measurements,
    const Acts::GeometryContext& /*gctx*/,
    Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::TrackStateProxy&
        trackState) const {
  const IndexSourceLink& sourceLink =
      trackState.getUncalibratedSourceLink().get<IndexSourceLink>();
  assert((sourceLink.index() < measurements.size()) and
         "Source link index is outside the container bounds");

  std::visit(
      [&trackState](const auto& meas) {
        trackState.allocateCalibrated(meas.size());
        trackState.setCalibrated(meas);
      },
      (measurements)[sourceLink.index()]);
}

ActsExamples::MeasurementCalibratorAdapter::MeasurementCalibratorAdapter(
    const MeasurementCalibrator& calibrator,
    const MeasurementContainer& measurements)
    : m_calibrator{calibrator}, m_measurements{measurements} {}

void ActsExamples::MeasurementCalibratorAdapter::calibrate(
    const Acts::GeometryContext& gctx,
    Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::TrackStateProxy
        trackState) const {
  return m_calibrator.calibrate(m_measurements, gctx, trackState);
}
