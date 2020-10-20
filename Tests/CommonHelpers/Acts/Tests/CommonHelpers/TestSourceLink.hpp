// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <cassert>
#include <iosfwd>

namespace Acts {
namespace Test {

/// A minimal source link implementation for testing.
///
/// Just stores a pointer to a measurement without actually linking back to a
/// detector readout. In addition it stores a source identifier that can be used
/// to store additional information. How this is interpreted depends on the
/// specific tests.
struct TestSourceLink {
  const FittableMeasurement<TestSourceLink>* measurement = nullptr;
  size_t sourceId = 0u;

  /// Construct from a measurement and optional source identifier.
  TestSourceLink(const FittableMeasurement<TestSourceLink>& m, size_t sid = 0u)
      : measurement(&m), sourceId(sid) {}
  /// Default-construct an invalid source link to satisfy SourceLinkConcept.
  TestSourceLink() = default;
  TestSourceLink(const TestSourceLink&) = default;
  TestSourceLink(TestSourceLink&&) = default;
  TestSourceLink& operator=(const TestSourceLink&) = default;
  TestSourceLink& operator=(TestSourceLink&&) = default;

  GeometryIdentifier geometryId() const {
    if (not measurement) {
      return GeometryIdentifier();
    }
    return MeasurementHelpers::getSurface(*measurement)->geometryId();
  }

  // implementing this a non-friend functions breaks the build. not sure why.
  friend bool operator==(const TestSourceLink& lhs, const TestSourceLink& rhs) {
    return (lhs.measurement == rhs.measurement) and
           (lhs.sourceId == rhs.sourceId);
  }
  friend bool operator!=(const TestSourceLink& lhs, const TestSourceLink& rhs) {
    return not(lhs == rhs);
  }
};

std::ostream& operator<<(std::ostream& os, const TestSourceLink& sourceLink);

/// Extract measurements from test source links.
///
/// The tests source links directly to a measurement so the calibration is
/// just a pass-through. Consequently, it does not depend on the track
/// parameters, but they still must be part of the interface.
struct TestSourceLinkCalibrator {
  /// Extract the measurement from a minimal source link.
  ///
  /// @tparam parameters_t Track parameters type
  /// @param sourceLink Input source link
  /// @param parameters Input track parameters (unused)
  template <typename parameters_t>
  const FittableMeasurement<TestSourceLink>& operator()(
      const TestSourceLink& sourceLink,
      const parameters_t& /* parameters */) const {
    assert(sourceLink.measurement and "Detected invalid minimal source link");
    return *sourceLink.measurement;
  }
};

}  // namespace Test
}  // namespace Acts
