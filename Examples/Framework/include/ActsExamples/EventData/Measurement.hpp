// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"

#include <cassert>
#include <vector>

namespace ActsExamples {

/// Variable measurement type that can contain all possible combinations.
using Measurement = ::Acts::BoundVariantMeasurement;
/// Container of measurements.
///
/// In contrast to the source links, the measurements themself must not be
/// orderable. The source links stored in the measurements are treated
/// as opaque here and no ordering is enforced on the stored measurements.
using MeasurementContainer = std::vector<Measurement>;

/// Calibrator to convert an index source link to a measurement.
class MeasurementCalibrator {
 public:
  /// Construct an invalid calibrator. Required to allow copying.
  MeasurementCalibrator() = default;
  /// Construct using a user-provided container to chose measurements from.
  MeasurementCalibrator(const MeasurementContainer& measurements)
      : m_measurements(&measurements) {}

  /// Find the measurement corresponding to the source link.
  ///
  /// @tparam parameters_t Track parameters type
  /// @param gctx The geometry context (unused)
  /// @param trackState The track state to calibrate
  void calibrate(
      const Acts::GeometryContext& /*gctx*/,
      Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::TrackStateProxy
          trackState) const {
    const IndexSourceLink& sourceLink =
        trackState.uncalibratedSourceLink().get<IndexSourceLink>();
    assert(m_measurements and
           "Undefined measurement container in DigitizedCalibrator");
    assert((sourceLink.index() < m_measurements->size()) and
           "Source link index is outside the container bounds");
    std::visit(
        [&trackState](const auto& meas) {
          trackState.allocateCalibrated(meas.size());
          trackState.setCalibrated(meas);
        },
        (*m_measurements)[sourceLink.index()]);
  }

 private:
  // use pointer so the calibrator is copyable and default constructible.
  const MeasurementContainer* m_measurements = nullptr;
};

}  // namespace ActsExamples
