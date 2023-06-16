// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <stdexcept>

namespace Acts {
namespace Test {

/// A minimal source link implementation for testing.
///
/// Instead of storing a reference to a measurement or raw data, the measurement
/// data is stored inline directly in the source link. Only 1d or 2d
/// measurements are supported to limit the overhead. Additionaly, a source
/// identifier is stored that can be used to store additional information. How
/// this is interpreted depends on the specific tests.
struct TestSourceLink final {
  GeometryIdentifier m_geometryId{};
  size_t sourceId = 0u;
  // use eBoundSize to indicate unused indices
  std::array<BoundIndices, 2> indices = {eBoundSize, eBoundSize};
  Acts::ActsVector<2> parameters;
  Acts::ActsSymMatrix<2> covariance;

  /// Construct a source link for a 1d measurement.
  TestSourceLink(BoundIndices idx, ActsScalar val, ActsScalar var,
                 GeometryIdentifier gid = GeometryIdentifier(), size_t sid = 0u)
      : m_geometryId(gid),
        sourceId(sid),
        indices{idx, eBoundSize},
        parameters(val, 0),
        covariance(Acts::ActsVector<2>(var, 0).asDiagonal()) {}
  /// Construct a source link for a 2d measurement.
  TestSourceLink(BoundIndices idx0, BoundIndices idx1,
                 const Acts::ActsVector<2>& params,
                 const Acts::ActsSymMatrix<2>& cov,
                 GeometryIdentifier gid = GeometryIdentifier(), size_t sid = 0u)
      : m_geometryId(gid),
        sourceId(sid),
        indices{idx0, idx1},
        parameters(params),
        covariance(cov) {}
  /// Default-construct an invalid source link to satisfy SourceLinkConcept.
  TestSourceLink() = default;
  TestSourceLink(const TestSourceLink&) = default;
  TestSourceLink(TestSourceLink&&) = default;
  TestSourceLink& operator=(const TestSourceLink&) = default;
  TestSourceLink& operator=(TestSourceLink&&) = default;

  constexpr size_t index() const { return sourceId; }

  GeometryIdentifier geometryId() const { return m_geometryId; }
};

bool operator==(const TestSourceLink& lhs, const TestSourceLink& rhs);
bool operator!=(const TestSourceLink& lhs, const TestSourceLink& rhs);
std::ostream& operator<<(std::ostream& os, const TestSourceLink& sourceLink);

/// Extract the measurement from a TestSourceLink.
///
/// @param gctx Unused
/// @param trackState TrackState to calibrated
/// @return The measurement used
template <typename trajectory_t>
Acts::BoundVariantMeasurement testSourceLinkCalibratorReturn(
    const GeometryContext& /*gctx*/,
    typename trajectory_t::TrackStateProxy trackState) {
  TestSourceLink sl =
      trackState.getUncalibratedSourceLink().template get<TestSourceLink>();
  if ((sl.indices[0] != Acts::eBoundSize) and
      (sl.indices[1] != Acts::eBoundSize)) {
    auto meas =
        makeMeasurement(trackState.getUncalibratedSourceLink(), sl.parameters,
                        sl.covariance, sl.indices[0], sl.indices[1]);
    trackState.allocateCalibrated(2);
    trackState.setCalibrated(meas);
    return meas;
  } else if (sl.indices[0] != Acts::eBoundSize) {
    auto meas = makeMeasurement(
        trackState.getUncalibratedSourceLink(), sl.parameters.head<1>(),
        sl.covariance.topLeftCorner<1, 1>(), sl.indices[0]);
    trackState.allocateCalibrated(1);
    trackState.setCalibrated(meas);
    return meas;
  } else {
    throw std::runtime_error(
        "Tried to extract measurement from invalid TestSourceLink");
  }
}
/// Extract the measurement from a TestSourceLink.
///
/// @param gctx Unused
/// @param trackState TrackState to calibrated
template <typename trajectory_t>
void testSourceLinkCalibrator(
    const GeometryContext& gctx,
    typename trajectory_t::TrackStateProxy trackState) {
  testSourceLinkCalibratorReturn<trajectory_t>(gctx, trackState);
}

}  // namespace Test
}  // namespace Acts
