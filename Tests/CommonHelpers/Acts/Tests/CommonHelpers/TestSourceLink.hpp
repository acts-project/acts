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
#include "Acts/EventData/SourceLinkConcept.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <cassert>
#include <iosfwd>

namespace Acts {

/// A minimal source link implementation.
///
/// Just stores a pointer to a measurement without actually linking back to a
/// detector readout. This is mostly intended as a simple test and demonstrator
/// type. Users are discouraged from using this. They are instead expected to
/// provide a custom source link implementation that is specific to their
/// requirements.
struct MinimalSourceLink {
  const FittableMeasurement<MinimalSourceLink>* measurement = nullptr;

  GeometryIdentifier geometryId() const {
    // invalid measurement results in an undefined geometry identifier
    if (not measurement) {
      return GeometryIdentifier();
    }
    return MeasurementHelpers::getSurface(*measurement)->geometryId();
  }
};

/// Compare two minimal source links for equality.
constexpr bool operator==(const MinimalSourceLink& lhs,
                          const MinimalSourceLink& rhs) {
  return lhs.measurement == rhs.measurement;
}

/// Compare two minimal source links for in-equality.
constexpr bool operator!=(const MinimalSourceLink& lhs,
                          const MinimalSourceLink& rhs) {
  return lhs.measurement != rhs.measurement;
}

std::ostream& operator<<(std::ostream& os, const MinimalSourceLink& sourceLink);

static_assert(SourceLinkConcept<MinimalSourceLink>,
              "MinimalSourceLink does not fulfill SourceLinkConcept");

/// Compute measurements from minimal source links.
///
/// The minimal source links directly to a measurement so the calibration is
/// just a pass-through. Consequently, it does not depend on the track
/// parameters, but they still must be part of the interface.
struct MinimalSourceLinkCalibrator {
  /// Compute the measurement from a minimal source link.
  ///
  /// @tparam parameters_t Track parameters type
  /// @param sourceLink Input source link
  /// @param parameters Input track parameters (unused)
  template <typename parameters_t>
  const FittableMeasurement<MinimalSourceLink>& operator()(
      const MinimalSourceLink& sourceLink,
      const parameters_t& /* parameters */) const {
    assert(sourceLink.measurement and "Detected invalid minimal source link");
    return *sourceLink.measurement;
  }
};

}  // namespace Acts
