// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/SourceLink.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include <ActsExamples/EventData/MeasurementCalibration.hpp>

#include <cassert>
#include <variant>

namespace Acts {
class VectorMultiTrajectory;
}  // namespace Acts

void ActsExamples::PassThroughCalibrator::calibrate(
    const MeasurementContainer& measurements,
    const ClusterContainer* /*clusters*/, const Acts::GeometryContext& /*gctx*/,
    const Acts::CalibrationContext& /*cctx*/,
    const Acts::SourceLink& sourceLink,
    Acts::VectorMultiTrajectory::TrackStateProxy& trackState) const {
  trackState.setUncalibratedSourceLink(sourceLink);
  const IndexSourceLink& idxSourceLink = sourceLink.get<IndexSourceLink>();

  assert((idxSourceLink.index() < measurements.size()) &&
         "Source link index is outside the container bounds");

  std::visit(
      [&trackState](const auto& meas) { trackState.setCalibrated(meas); },
      measurements[idxSourceLink.index()]);
}

ActsExamples::MeasurementCalibratorAdapter::MeasurementCalibratorAdapter(
    const MeasurementCalibrator& calibrator,
    const MeasurementContainer& measurements, const ClusterContainer* clusters)
    : m_calibrator{calibrator},
      m_measurements{measurements},
      m_clusters{clusters} {}

void ActsExamples::MeasurementCalibratorAdapter::calibrate(
    const Acts::GeometryContext& gctx, const Acts::CalibrationContext& cctx,
    const Acts::SourceLink& sourceLink,
    Acts::VectorMultiTrajectory::TrackStateProxy trackState) const {
  return m_calibrator.calibrate(m_measurements, m_clusters, gctx, cctx,
                                sourceLink, trackState);
}
