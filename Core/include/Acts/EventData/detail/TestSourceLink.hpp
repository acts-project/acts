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
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <stdexcept>

namespace Acts::detail::Test {

/// A minimal source link implementation for testing.
///
/// Instead of storing a reference to a measurement or raw data, the measurement
/// data is stored inline directly in the source link. Only 1d or 2d
/// measurements are supported to limit the overhead. Additionally, a source
/// identifier is stored that can be used to store additional information. How
/// this is interpreted depends on the specific tests.
struct TestSourceLink final {
  GeometryIdentifier m_geometryId{};
  std::size_t sourceId = 0u;
  // use eBoundSize to indicate unused indices
  std::array<BoundIndices, 2> indices = {eBoundSize, eBoundSize};
  Acts::ActsVector<2> parameters;
  Acts::ActsSquareMatrix<2> covariance;

  /// Construct a source link for a 1d measurement.
  TestSourceLink(BoundIndices idx, ActsScalar val, ActsScalar var,
                 GeometryIdentifier gid = GeometryIdentifier(),
                 std::size_t sid = 0u)
      : m_geometryId(gid),
        sourceId(sid),
        indices{idx, eBoundSize},
        parameters(val, 0),
        covariance(Acts::ActsVector<2>(var, 0).asDiagonal()) {}
  /// Construct a source link for a 2d measurement.
  TestSourceLink(BoundIndices idx0, BoundIndices idx1,
                 const Acts::ActsVector<2>& params,
                 const Acts::ActsSquareMatrix<2>& cov,
                 GeometryIdentifier gid = GeometryIdentifier(),
                 std::size_t sid = 0u)
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
  bool operator==(const TestSourceLink& rhs) const {
    return (m_geometryId == rhs.m_geometryId) && (sourceId == rhs.sourceId) &&
           (indices == rhs.indices) && (parameters == rhs.parameters) &&
           (covariance == rhs.covariance);
  }
  bool operator!=(const TestSourceLink& rhs) const { return !(*this == rhs); }
  std::ostream& print(std::ostream& os) const {
    os << "TestsSourceLink(geometryId=" << m_geometryId
       << ",sourceId=" << sourceId;
    if (indices[0] != eBoundSize) {
      os << ",index0=" << indices[0];
    }
    if (indices[1] != eBoundSize) {
      os << ",index1=" << indices[1];
    }
    os << ")";
    return os;
  }
  constexpr std::size_t index() const { return sourceId; }

  struct SurfaceAccessor {
    const Acts::TrackingGeometry& trackingGeometry;

    const Acts::Surface* operator()(const Acts::SourceLink& sourceLink) const {
      const auto& testSourceLink = sourceLink.get<TestSourceLink>();
      return trackingGeometry.findSurface(testSourceLink.m_geometryId);
    }
  };
};

inline std::ostream& operator<<(std::ostream& os,
                                const TestSourceLink& sourceLink) {
  return sourceLink.print(os);
}

/// Extract the measurement from a TestSourceLink.
///
/// @param gctx Unused
/// @param trackState TrackState to calibrated
/// @return The measurement used
template <typename trajectory_t>
Acts::BoundVariantMeasurement testSourceLinkCalibratorReturn(
    const GeometryContext& /*gctx*/, const CalibrationContext& /*cctx*/,
    const SourceLink& sourceLink,
    typename trajectory_t::TrackStateProxy trackState) {
  TestSourceLink sl = sourceLink.template get<TestSourceLink>();

  trackState.setUncalibratedSourceLink(sourceLink);

  if ((sl.indices[0] != Acts::eBoundSize) &&
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
    const GeometryContext& gctx, const CalibrationContext& cctx,
    const SourceLink& sourceLink,
    typename trajectory_t::TrackStateProxy trackState) {
  testSourceLinkCalibratorReturn<trajectory_t>(gctx, cctx, sourceLink,
                                               trackState);
}

}  // namespace Acts::detail::Test
